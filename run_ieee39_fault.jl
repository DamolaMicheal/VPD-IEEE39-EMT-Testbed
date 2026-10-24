#!/usr/bin/env julia
#IEEE 39-bus three-phase fault study.
#Run this file to solve the EMT model and export voltage/angle data.
using PowerSystems
using PowerSimulationsDynamics
using PowerFlows
using Sundials
using Plots
using LinearAlgebra
using Statistics
using Dates
using Printf
using Logging
using CSV
using DataFrames
using SparseArrays
const PSY  = PowerSystems
const PSID = PowerSimulationsDynamics

include(joinpath(@__DIR__, "case39_matpower.jl"))
include(joinpath(@__DIR__, "emt_parameters.jl"))

#Load-change events for the selected scenario.  Each event gives the
#reference value of the perturbation load at that instant; the realized
#simulation steps are the differences between consecutive references.
function scenario_load_events()
    sp = SCENARIO
    events = Vector{Dict{Symbol,Any}}()
    add_event!(bus, t, p, q) = push!(events, Dict(:buses => [bus], :event_time_s => t,
                                                   :p_pu_map => Dict(bus => p), :q_pu_map => Dict(bus => q)))
    base_p = Dict{Int,Float64}()
    if sp.p_step != 0.0
        add_event!(sp.bus_step, sp.t_step, sp.p_step, LOAD_STEP_Q_PU)
        base_p[sp.bus_step] = sp.p_step
    end
    if sp.window !== nothing
        t0, t1 = sp.window
        base = get(base_p, sp.bus_osc, 0.0)
        n = round(Int, (t1 - t0) / SINE_UPDATE_DT)
        abs(t0 + n * SINE_UPDATE_DT - t1) < 1e-9 || error("SINE_UPDATE_DT must divide the oscillation window")
        for k in 0:n-1
            tk = t0 + k * SINE_UPDATE_DT
            add_event!(sp.bus_osc, tk, base + sp.amp * sin(2π * sp.f_hz * (tk - t0)), LOAD_STEP_Q_PU + SINE_Q_PU)
        end
        add_event!(sp.bus_osc, t1, base, LOAD_STEP_Q_PU)   # oscillation off, underlying load unchanged
    end
    return events
end

function describe_scenario()
    sp = SCENARIO
    parts = String[]
    sp.p_step != 0.0 && push!(parts, "$(sp.p_step) pu P step at bus $(sp.bus_step) at $(sp.t_step) s (kept on)")
    if sp.window !== nothing
        t0, t1 = sp.window
        push!(parts, "$(sp.amp) pu, $(sp.f_hz) Hz P oscillation at bus $(sp.bus_osc) from $t0 s to $t1 s")
    end
    if BUS_FAULT_ENABLE
        if FAULT_SCENARIO == :line_fault_line_trip
            push!(parts, "3-phase-to-ground fault at $(100 * FAULT_LOCATION_FROM)% of line $FAULT_LINE_FROM-$FAULT_LINE_TO " *
                         "at $FAULT_APPLY_TIME_S s, cleared at $FAULT_CLEAR_TIME_S s by tripping line $FAULT_LINE_FROM-$FAULT_LINE_TO")
        elseif FAULT_SCENARIO == :bus_fault_line_trip
            push!(parts, "3-phase fault at bus $BUS_FAULT_BUS at $FAULT_APPLY_TIME_S s, " *
                         "cleared at $BUS_FAULT_CLEAR_TIME_S s by tripping line $BUS_FAULT_TRIP_LINE_FROM-$BUS_FAULT_TRIP_LINE_TO")
        end
    end
    isempty(parts) && push!(parts, "no disturbance")
    return join(parts, "; ")
end

# Asset and parameter helpers
is_sg_bus(bus::Int)  = bus in SG_BUSES
is_gfm_bus(bus::Int) = bus in GFM_BUSES
is_gfl_bus(bus::Int) = bus in GFL_BUSES

function model_label(bus::Int)
    if is_sg_bus(bus)
        return "SG"
    elseif is_gfm_bus(bus)
        return "GFM"
    elseif is_gfl_bus(bus)
        return "GFL"
    else
        return "Unknown"
    end
end

function getH(bus::Int)
    haskey(H_ALL_map_ext, bus) && return H_ALL_map_ext[bus]
    @warn "No inertia H found for bus; using default H = 4.0" bus = bus model = model_label(bus)
    return 4.0
end

function getDsg(bus::Int)
    haskey(D_SG_map_ext, bus) && return D_SG_map_ext[bus]
    @warn "No SG damping D found for bus; using D_SG_DEFAULT" bus = bus model = model_label(bus) D_default = D_SG_DEFAULT
    return D_SG_DEFAULT
end

function getGfmTa(bus::Int)
    haskey(GFM_TA_map_ext, bus) && return GFM_TA_map_ext[bus]
    @warn "No GFM Ta found for bus; using GFM_TA_DEFAULT" bus = bus model = model_label(bus) Ta_default = GFM_TA_DEFAULT
    return GFM_TA_DEFAULT
end

function getGfmKw(bus::Int)
    haskey(GFM_KW_map_ext, bus) && return GFM_KW_map_ext[bus]
    @warn "No GFM Kω found for bus; using GFM_KW_DEFAULT" bus = bus model = model_label(bus) Kω_default = GFM_KW_DEFAULT
    return GFM_KW_DEFAULT
end

function getGfmKd(bus::Int)
    haskey(GFM_KD_map_ext, bus) && return GFM_KD_map_ext[bus]
    @warn "No GFM kd found for bus; using GFM_KD_DEFAULT" bus = bus model = model_label(bus) kd_default = GFM_KD_DEFAULT
    return GFM_KD_DEFAULT
end

function getGflKw(bus::Int)
    haskey(GFL_KW_map_ext, bus) && return GFL_KW_map_ext[bus]
    @warn "No GFL Kω found for bus; using GFL_KW_DEFAULT" bus = bus model = model_label(bus) Kω_default = GFL_KW_DEFAULT
    return GFL_KW_DEFAULT
end

#Validation
function validate_partitions!(sg::Vector{Int}, gfm::Vector{Int}, gfl::Vector{Int})
    sgset  = Set(sg)
    gfmset = Set(gfm)
    gflset = Set(gfl)

    if 39 ∉ sgset
        error("Bus 39 must be in SG_BUSES.")
    end
    if (39 ∈ gfmset) || (39 ∈ gflset)
        error("Bus 39 must NOT be in GFM_BUSES or GFL_BUSES.")
    end

    ov1 = intersect(sgset, gfmset)
    ov2 = intersect(sgset, gflset)
    ov3 = intersect(gfmset, gflset)
    if !isempty(ov1) || !isempty(ov2) || !isempty(ov3)
        error("Partitions overlap. SG∩GFM=$(collect(ov1)), SG∩GFL=$(collect(ov2)), GFM∩GFL=$(collect(ov3))")
    end

    allowed = Set(30:39)
    if !issubset(sgset, allowed) || !issubset(gfmset, allowed) || !issubset(gflset, allowed)
        error("SG_BUSES/GFM_BUSES/GFL_BUSES must be within 30..39.")
    end

    allowed_gfl = Set(30:38)
    if !issubset(gflset, allowed_gfl)
        bad = collect(setdiff(gflset, allowed_gfl))
        error("GFL_BUSES must be within 30..38 only. Bad=$(bad)")
    end

    cover = union(union(sgset, gfmset), gflset)
    if cover != allowed
        missing = collect(setdiff(allowed, cover))
        extra   = collect(setdiff(cover, allowed))
        error("Partitions must cover ALL buses 30..39 exactly. Missing=$(missing), Extra=$(extra)")
    end

    return true
end

# Load-event helpers
struct LoadEventSpec
    buses::Vector{Int}
    event_time_s::Float64
    p_pu_map::Dict{Int,Float64}
    q_pu_map::Dict{Int,Float64}
end

struct RealizedLoadStep
    bus::Int
    t_step::Float64
    ΔP_MW::Float64
    ΔQ_MVAr::Float64
    refP_MW::Float64
    refQ_MVAr::Float64
    refP_pu::Float64
    refQ_pu::Float64
    deltaP_pu::Float64
    deltaQ_pu::Float64
end

function _event_p_map(ev)
    if haskey(ev, :p_pu_map)
        return ev[:p_pu_map]
    elseif haskey(ev, :pu_map)
        return ev[:pu_map]
    end
    return Dict{Int,Float64}()
end

function _event_q_map(ev)
    haskey(ev, :q_pu_map) && return ev[:q_pu_map]
    return Dict{Int,Float64}()
end

function _to_int_float_dict(raw_map)
    out = Dict{Int,Float64}()
    for (k, v) in raw_map
        out[Int(k)] = Float64(v)
    end
    return out
end

function load_event_specs_from_config()
    specs = LoadEventSpec[]
    for ev in scenario_load_events()
        buses = Vector{Int}(ev[:buses])
        t     = Float64(ev[:event_time_s])
        p_map = _to_int_float_dict(_event_p_map(ev))
        q_map = _to_int_float_dict(_event_q_map(ev))
        push!(specs, LoadEventSpec(buses, t, p_map, q_map))
    end
    sort!(specs, by = x -> x.event_time_s)
    return specs
end

function build_realized_load_steps_from_events()
    specs = load_event_specs_from_config()
    prevP_MW   = Dict{Int,Float64}()
    prevQ_MVAr = Dict{Int,Float64}()
    realized_steps = RealizedLoadStep[]
    for ev in specs
        for b in ev.buses
            refP_pu = get(ev.p_pu_map, b, 0.0)
            refQ_pu = get(ev.q_pu_map, b, 0.0)
            refP_MW   = refP_pu * SYS_BASE_MVA
            refQ_MVAr = refQ_pu * SYS_BASE_MVA
            oldP_MW   = get(prevP_MW, b, 0.0)
            oldQ_MVAr = get(prevQ_MVAr, b, 0.0)
            ΔP_MW   = refP_MW - oldP_MW
            ΔQ_MVAr = refQ_MVAr - oldQ_MVAr
            prevP_MW[b]   = refP_MW
            prevQ_MVAr[b] = refQ_MVAr
            push!(realized_steps,
                RealizedLoadStep(b, ev.event_time_s, ΔP_MW, ΔQ_MVAr, refP_MW, refQ_MVAr,
                                 refP_pu, refQ_pu, ΔP_MW / SYS_BASE_MVA, ΔQ_MVAr / SYS_BASE_MVA))
        end
    end
    sort!(realized_steps, by = x -> x.t_step)
    return realized_steps
end

function print_load_change_schedule(realized_steps::Vector{RealizedLoadStep})
    println("\nDisturbance scenario: $DISTURBANCE_SCENARIO")
    println("  " * describe_scenario())
    println("  realized load-change events: $(length(realized_steps))")
    if length(realized_steps) <= 12
        for st in realized_steps
            @printf("  t = %.4f s  bus %d  P* = %+.4f pu (step %+.4f)  Q* = %+.4f pu\n",
                    st.t_step, st.bus, st.refP_pu, st.deltaP_pu, st.refQ_pu)
        end
    else
        st = realized_steps[1]; en = realized_steps[end]
        @printf("  first: t = %.4f s  P* = %+.4f pu;  last: t = %.4f s  P* = %+.4f pu\n",
                st.t_step, st.refP_pu, en.t_step, en.refP_pu)
    end
end

event_times_from_realized_steps(realized_steps::Vector{RealizedLoadStep}) =
    unique(sort([s.t_step for s in realized_steps]))

# Dynamic line setup
# Replace each transformer by a unit-ratio RL Line so it can become a DynamicBranch.
function transformers_to_dynamic_lines!(sys::System)
    TRANSFORMERS_DYNAMIC || return 0
    n = 0
    for T in (Transformer2W, TapTransformer)
        for xf in collect(get_components(T, sys))
            arc = PSY.get_arc(xf)
            rating = try Float64(PSY.get_rating(xf)) catch; 10.0 end
            line = Line(name = "dynxf_" * PSY.get_name(xf), available = true,
                        active_power_flow = 0.0, reactive_power_flow = 0.0, arc = arc,
                        r = Float64(PSY.get_r(xf)), x = Float64(PSY.get_x(xf)),
                        b = (from = TRANSFORMER_DYN_B / 2, to = TRANSFORMER_DYN_B / 2),
                        rating = rating, angle_limits = (min = -pi / 2, max = pi / 2))
            remove_component!(sys, xf)
            # PSY refuses a Line between buses of different base voltage (the two sides
            # of a transformer); the per-unit model is unaffected, so skip that check.
            add_component!(sys, line; skip_validation = true)
            n += 1
        end
    end
    println("Replaced $n transformers by unit-ratio RL lines with b = $(TRANSFORMER_DYN_B) pu (TRANSFORMERS_DYNAMIC).")
    return n
end

function add_dynamic_lines!(sys::System)
    n_added = 0
    for line in get_components(Line, sys)
        try
            dyn_branch = DynamicBranch(line)
            add_component!(sys, dyn_branch)
            n_added += 1
        catch e
            @warn "Skipping line -> DynamicBranch conversion (constructor/add failed)" line = PSY.get_name(line) err = string(e)
        end
    end
    println("Converted $n_added AC lines to dynamic RL branches.")
    return n_added
end

function gen_at_bus(sys::System, busnum::Int)
    for g in get_components(Generator, sys)
        if PSY.get_number(PSY.get_bus(g)) == busnum
            return g
        end
    end
    error("No Generator found at bus $busnum")
end

function update_system_voltages!(sys::System, pf_data::PowerFlowData)
    for bus in get_components(Bus, sys)
        bus_num = PSY.get_number(bus)
        if haskey(pf_data.bus_lookup, bus_num)
            idx = pf_data.bus_lookup[bus_num]
            PSY.set_magnitude!(bus, pf_data.bus_magnitude[idx])
            PSY.set_angle!(bus, pf_data.bus_angles[idx])
        end
    end
end

function enforce_load_bases!(sys::System; base_mva::Float64 = SYS_BASE_MVA)
    for ld in get_components(ElectricLoad, sys)
        try
            PSY.set_base_power!(ld, base_mva)
        catch
        end
    end
end

#Give zero-resistance branches a small r before the power flow so every edge has
#r/x > 0.  The exported network model applies the same rule, so plant and
#estimator agree.  Tapped transformers are TapTransformer objects, so they are included.
function enforce_min_line_r!(sys::System; rmin::Float64 = LINE_RMIN, xr_cap::Float64 = LINE_XR_CAP)
    n = 0
    for component_type in (Line, Transformer2W, TapTransformer)
        for branch in get_components(component_type, sys)
            r = Float64(PSY.get_r(branch))
            x = Float64(PSY.get_x(branch))
            target = if r > 1e-12
                ENABLE_LINE_XR_CAP ? max(r, abs(x) / xr_cap) : r
            else
                ENABLE_LINE_XR_CAP ? max(rmin, abs(x) / xr_cap) : rmin
            end
            if target > r + 1e-12
                # Do not silently continue with inconsistent PF/VPD resistance.
                PSY.set_r!(branch, target)
                n += 1
            end
        end
    end
    println("Series resistance regularized on $n Lines/Transformers (rmin=$rmin, X/R cap=$(ENABLE_LINE_XR_CAP)).")
    return n
end

# Use passive impedance loads except perturbation loads
function loads_to_constant_impedance!(sys::System)
    converted = 0
    for ld in collect(get_components(PowerLoad, sys))
        nm = PSY.get_name(ld)
        startswith(nm, "PerturbLoad") && continue
        b = PSY.get_bus(ld)
        V = max(try PSY.get_magnitude(b) catch; 1.0 end, 0.5)
        P = try PSY.get_active_power(ld)   catch; 0.0 end
        Q = try PSY.get_reactive_power(ld) catch; 0.0 end
        (abs(P) < 1e-9 && abs(Q) < 1e-9) && continue
        Y = complex(P, -Q) / V^2
        adm = nothing
        try
            adm = FixedAdmittance(name = "Zload_" * nm, available = true, bus = b, Y = Y)
        catch
            adm = FixedAdmittance("Zload_" * nm, true, b, Y)
        end
        remove_component!(sys, ld)
        add_component!(sys, adm)
        converted += 1
    end
    println("Converted $converted constant-power loads to constant impedance.")
    return converted
end

function sync_static_gens_to_pf!(sys::System, pf_data::PowerFlowData; buses::Vector{Int})
    for b in buses
        if !haskey(pf_data.bus_lookup, b)
            @warn "Bus not found in PF bus_lookup; skipping" bus = b
            continue
        end
        idx = pf_data.bus_lookup[b]
        Pinj = try pf_data.bus_activepower_injection[idx] catch; 0.0 end
        Qinj = try pf_data.bus_reactivepower_injection[idx] catch; 0.0 end
        Pwd  = try pf_data.bus_activepower_withdrawals[idx] catch; 0.0 end
        Qwd  = try pf_data.bus_reactivepower_withdrawals[idx] catch; 0.0 end
        Pnet_pu = Pinj - Pwd
        Qnet_pu = Qinj - Qwd
        g = gen_at_bus(sys, b)
        PSY.set_active_power!(g, Pnet_pu)
        try
            PSY.set_reactive_power!(g, Qnet_pu)
        catch
        end
    end
end

#Controller reference helpers
function get_generator_pref_qref(g::Generator)
    p_ref = try
        PSY.get_active_power(g)
    catch
        0.0
    end
    q_ref = try
        PSY.get_reactive_power(g)
    catch
        0.0
    end

#Convert controller references to device base
    mbase = try PSY.get_base_power(g) catch; SYS_BASE_MVA end
    scale = SYS_BASE_MVA / mbase
    return p_ref * scale, q_ref * scale
end

function get_bus_vref(sys::System, g::Generator)
    bus = PSY.get_bus(g)
    v_ref = try
        PSY.get_magnitude(bus)
    catch
        1.0
    end
    return v_ref
end

function print_inverter_reference_summary(sys::System; gfm_buses::Vector{Int}, gfl_buses::Vector{Int})
    println("\n" * "="^100)
    println("INVERTER OUTER-CONTROL REFERENCE VALUES (from initialized PF operating point)")
    println("="^100)
    if !isempty(gfm_buses)
        println("GFM buses:")
        println(rpad("Bus", 8), rpad("P_ref [pu]", 18), rpad("Q_ref [pu]", 18), rpad("V_ref [pu]", 18))
        println("-"^100)
        for b in gfm_buses
            g = gen_at_bus(sys, b)
            p_ref0, q_ref0 = get_generator_pref_qref(g)
            v_ref0 = get_bus_vref(sys, g)
            println(rpad(string(b), 8), rpad(@sprintf("%.6f", p_ref0), 18),
                    rpad(@sprintf("%.6f", q_ref0), 18), rpad(@sprintf("%.6f", v_ref0), 18))
        end
        println("-"^100)
    end
    if !isempty(gfl_buses)
        println("GFL buses:")
        println(rpad("Bus", 8), rpad("P_ref [pu]", 18), rpad("Q_ref [pu]", 18))
        println("-"^100)
        for b in gfl_buses
            g = gen_at_bus(sys, b)
            p_ref0, q_ref0 = get_generator_pref_qref(g)
            println(rpad(string(b), 8), rpad(@sprintf("%.6f", p_ref0), 18), rpad(@sprintf("%.6f", q_ref0), 18))
        end
        println("-"^100)
    end
    println("="^100 * "\n")
end

function print_dynamic_reference_summary(sys::System; buses::Vector{Int})
    println("\n" * "="^112)
    println("DYNAMIC MODEL REFERENCE SUMMARY AFTER PF SYNC AND mBase SETTING")
    println("="^112)
    @printf("%-8s %-10s %-16s %-16s %-16s %-14s %-14s\n",
            "Bus", "Model", "P_ref[pu]", "Q_ref[pu]", "V_ref[pu]", "mBase[MVA]", "P_ref[MW]")
    println("-"^112)
    for b in buses
        g = gen_at_bus(sys, b)
        p_ref0, q_ref0 = get_generator_pref_qref(g)
        v_ref0 = get_bus_vref(sys, g)
        mbase = try PSY.get_base_power(g) catch; NaN end
        @printf("%-8d %-10s %-16.6f %-16.6f %-16.6f %-14.2f %-14.3f\n",
                b, model_label(b), p_ref0, q_ref0, v_ref0, mbase, p_ref0 * mbase)
    end
    println("="^112 * "\n")
end

# Power-flow export
function _pf_angle_to_deg(ang::Real)
    if abs(ang) <= 2π + 0.5
        return ang * (180.0 / π)
    else
        return ang
    end
end

function build_pf_bus_dataframe(pf_data::PowerFlowData; sys_base_mva::Float64 = SYS_BASE_MVA)
    busnums = sort(collect(keys(pf_data.bus_lookup)))
    rows = Vector{NamedTuple}(undef, length(busnums))
    for (k, b) in enumerate(busnums)
        idx = pf_data.bus_lookup[b]
        Vm = try pf_data.bus_magnitude[idx] catch; NaN end
        Va = try pf_data.bus_angles[idx]    catch; NaN end
        Pinj = try pf_data.bus_activepower_injection[idx]   catch; 0.0 end
        Qinj = try pf_data.bus_reactivepower_injection[idx] catch; 0.0 end
        Pwd  = try pf_data.bus_activepower_withdrawals[idx] catch; 0.0 end
        Qwd  = try pf_data.bus_reactivepower_withdrawals[idx] catch; 0.0 end
        Pnet = Pinj - Pwd
        Qnet = Qinj - Qwd
        rows[k] = (
            bus = b, Vm_pu = Vm, Va_deg = _pf_angle_to_deg(Va),
            Pinj_pu = Pinj, Qinj_pu = Qinj, Pwd_pu = Pwd, Qwd_pu = Qwd,
            Pnet_pu = Pnet, Qnet_pu = Qnet,
            Pinj_MW = Pinj * sys_base_mva, Qinj_MVAr = Qinj * sys_base_mva,
            Pwd_MW = Pwd * sys_base_mva, Qwd_MVAr = Qwd * sys_base_mva,
            Pnet_MW = Pnet * sys_base_mva, Qnet_MVAr = Qnet * sys_base_mva
        )
    end
    return DataFrame(rows)
end


#Operating-point table for Python: the solved power flow plus the EMT's own
#pre-disturbance phasor at every bus (Vm_emt0_pu, Va_emt0_deg).  The EMT settles a
#few 1e-4 pu away from the power flow at load buses, and the VPD has to be
#linearized at, and plotted on, the EMT level.
function operating_point_table(pf_data::PowerFlowData, net)
    df = build_pf_bus_dataframe(pf_data)
    idx = Dict(b => i for (i, b) in enumerate(net.buses))
    V0 = [haskey(idx, b) ? net.V0[idx[b]] : complex(NaN, NaN) for b in df.bus]
    df.Vm_emt0_pu = abs.(V0)
    df.Va_emt0_deg = rad2deg.(angle.(V0))
    dv = abs.(df.Vm_emt0_pu .- df.Vm_pu)
    k = argmax(replace(dv, NaN => -Inf))
    @printf("Operating point: EMT pre-disturbance |V| differs from the power flow by up to %.2e pu (bus %d)\n",
            dv[k], df.bus[k])
    return df
end

#Series impedances of every branch as seen by the simulated system, keyed by
#the normalized (from, to) pair.  Used only to cross-check the exported model.
function system_branch_impedances(sys::System)
    out = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for T in (Line, DynamicBranch, Transformer2W, TapTransformer)
        comps = try collect(get_components(T, sys)) catch; [] end
        for c in comps
            inner = try PSY.get_branch(c) catch; c end
            pair = try
                line_pair(inner)
            catch
                try
                    ft = vpd_dynamic_branch_endpoints(c)
                    normalize_pair(ft[1], ft[2])
                catch
                    nothing
                end
            end
            pair === nothing && continue
            r = _line_parameter_or_default(inner, :get_r, NaN)
            x = _line_parameter_or_default(inner, :get_x, NaN)
            (isfinite(r) && isfinite(x)) || continue
            out[pair] = (r, x)
        end
    end
    return out
end

#Network file for Python, with the same r floor as the simulated system; every
#row is checked against the components actually in the system.
function write_vpd_network_model_csv(sys::System, outdir::String, run_id::String)
    sys_rx = system_branch_impedances(sys)
    rows = NamedTuple[]
    n_checked = 0
    for (f, t, r, x, b, tap_raw) in VPD_BRANCH_TABLE
        effective_r = if r > 1e-12
            ENABLE_LINE_XR_CAP ? max(r, abs(x) / LINE_XR_CAP) : r
        else
            ENABLE_LINE_XR_CAP ? max(LINE_RMIN, abs(x) / LINE_XR_CAP) : LINE_RMIN
        end
        pair = normalize_pair(f, t)
        if haskey(sys_rx, pair)
            rs, xs = sys_rx[pair]
            (abs(rs - effective_r) <= 1e-9 && abs(xs - x) <= 1e-9) ||
                error("Exported network model differs from the simulated system on branch $pair: " *
                      "export r=$effective_r x=$x, system r=$rs x=$xs")
            n_checked += 1
        end
        is_transformer = abs(tap_raw) > 1e-12
        push!(rows, (
            run_id = run_id, from_bus = f, to_bus = t,
            r_pu = effective_r, x_pu = x, b_pu = b,
            tap = is_transformer ? tap_raw : 0.0, shift_deg = 0.0, status = 1,
            branch_kind = is_transformer ? "algebraic_transformer" : "dynamic_rl_line",
        ))
    end
    count(r -> r.branch_kind == "dynamic_rl_line", rows) == VPD_EXPECTED_DYNAMIC_BRANCHES ||
        error("Expected $(VPD_EXPECTED_DYNAMIC_BRANCHES) dynamic lines")
    count(r -> r.branch_kind == "algebraic_transformer", rows) == VPD_EXPECTED_STATIC_TRANSFORMERS ||
        error("Expected $(VPD_EXPECTED_STATIC_TRANSFORMERS) algebraic transformers")
    n_checked > 0 || @warn "Could not read branch impedances from the system; network model exported unchecked"
    @info "VPD network model exported" branches = length(rows) checked_against_system = n_checked
    name = @sprintf("vpd_network_model_%s.csv", run_id)
    CSV.write(joinpath(outdir, name), DataFrame(rows))
    return name
end

function bus_component(sys::System, busnum::Int)
    for b in get_components(Bus, sys)
        if PSY.get_number(b) == busnum
            return b
        end
    end
    error("Bus $busnum not found in system")
end

function ensure_perturbable_load_on_bus!(sys::System, busnum::Int; base_mva::Float64 = SYS_BASE_MVA)
    name = @sprintf("PerturbLoad_bus%02d", busnum)
    for ld in get_components(ElectricLoad, sys)
        if PSY.get_name(ld) == name
            try PSY.set_base_power!(ld, base_mva) catch end
            return ld
        end
    end
    b = bus_component(sys, busnum)
    ld = nothing
    try
        ld = PowerLoad(name = name, available = true, bus = b, active_power = 0.0,
                       reactive_power = 0.0, base_power = base_mva,
                       max_active_power = 10.0, max_reactive_power = 10.0)
    catch
        try
            ld = PowerLoad(name, true, b, 0.0, 0.0, base_mva, 10.0, 10.0)
        catch
            ld = PowerLoad(name, true, b, 0.0, 0.0, base_mva)
        end
    end
    add_component!(sys, ld)
    try PSY.set_base_power!(ld, base_mva) catch end
    return ld
end

function ensure_perturbable_loads_all_buses!(sys::System; base_mva::Float64 = SYS_BASE_MVA)
    for busnum in 1:39
        ensure_perturbable_load_on_bus!(sys, busnum; base_mva = base_mva)
    end
end

function perturb_load_at_bus(sys::System, busnum::Int)
    target = @sprintf("PerturbLoad_bus%02d", busnum)
    for ld in get_components(ElectricLoad, sys)
        if PSY.get_name(ld) == target
            if PSY.get_number(PSY.get_bus(ld)) != busnum
                error("Found $target but attached to bus $(PSY.get_number(PSY.get_bus(ld))) not $busnum")
            end
            return ld
        end
    end
    error("Perturb load $target not found")
end

function print_loads_on_bus(sys::System, busnum::Int)
    println("\nLoads on bus $busnum:")
    found = false
    for ld in get_components(ElectricLoad, sys)
        if PSY.get_number(PSY.get_bus(ld)) == busnum
            found = true
            basep = try PSY.get_base_power(ld) catch; NaN end
            ppu   = try PSY.get_active_power(ld) catch; NaN end
            qpu   = try PSY.get_reactive_power(ld) catch; NaN end
            println("  - $(PSY.get_name(ld)) :: $(typeof(ld)) | base=$(basep) | P=$(ppu) pu | Q=$(qpu) pu")
        end
    end
    !found && println("  (none)")
end

function make_load_change_constantP(sys::System; bus::Int, t_step::Float64, ΔP_MW::Float64)
    ld = perturb_load_at_bus(sys, bus)
    ΔP_pu = ΔP_MW / SYS_BASE_MVA
    for field in (:active_power, :P, :p, :P_ref_power, :P_ref)
        try
            return LoadChange(t_step, ld, field, ΔP_pu)
        catch
        end
    end
    error("Could not apply active-power load step: no recognized active-power field found on $(typeof(ld)).")
end

function make_load_change_constantQ(sys::System; bus::Int, t_step::Float64, ΔQ_MVAr::Float64)
    ld = perturb_load_at_bus(sys, bus)
    ΔQ_pu = ΔQ_MVAr / SYS_BASE_MVA
    for field in (:reactive_power, :Q, :q, :Q_ref_power, :Q_ref)
        try
            return LoadChange(t_step, ld, field, ΔQ_pu)
        catch
        end
    end
    error("Could not apply reactive-power load step: no recognized reactive-power field found on $(typeof(ld)).")
end

function make_load_changes_from_event_schedule(sys::System, realized_steps::Vector{RealizedLoadStep})
    perts = PowerSimulationsDynamics.Perturbation[]
    for s in realized_steps
        if abs(s.ΔP_MW) > 1e-12
            push!(perts, make_load_change_constantP(sys; bus = s.bus, t_step = s.t_step, ΔP_MW = s.ΔP_MW))
        end
        if abs(s.ΔQ_MVAr) > 1e-12
            push!(perts, make_load_change_constantQ(sys; bus = s.bus, t_step = s.t_step, ΔQ_MVAr = s.ΔQ_MVAr))
        end
    end
    return perts
end

# Network-switch fault helpers
function print_bus_fault_line_trip_schedule()
    println("\n" * "="^118)
    println("LINE-FAULT + LINE-TRIP DISTURBANCE SCHEDULE")
    println("="^118)
    @printf("Fault-on   : t = %.3f s, 3φ-G fault on line %d--%d at %.1f%% from bus %d, Zf = %.6g %+.6gim pu\n",
            FAULT_APPLY_TIME_S, FAULT_LINE_FROM, FAULT_LINE_TO,
            100.0 * FAULT_LOCATION_FROM, FAULT_LINE_FROM, FAULT_R_PU, FAULT_X_PU)
    @printf("Fault-clear: t = %.3f s, remove fault and trip/open line %d--%d\n",
            FAULT_CLEAR_TIME_S, FAULT_LINE_FROM, FAULT_LINE_TO)
    @printf("Fault duration = %.3f s = %.1f ms\n",
            FAULT_CLEAR_TIME_S - FAULT_APPLY_TIME_S, 1000.0 * (FAULT_CLEAR_TIME_S - FAULT_APPLY_TIME_S))
    println("="^118 * "\n")
end


function selected_fault_clear_time()
    if FAULT_SCENARIO == :line_fault_line_trip
        return FAULT_CLEAR_TIME_S
    elseif FAULT_SCENARIO == :bus_fault_line_trip
        return BUS_FAULT_CLEAR_TIME_S
    else
        error("Unsupported FAULT_SCENARIO=$(FAULT_SCENARIO). Use :line_fault_line_trip or :bus_fault_line_trip.")
    end
end

function fault_event_times()
    BUS_FAULT_ENABLE || return Float64[]
    return unique(sort([FAULT_APPLY_TIME_S, selected_fault_clear_time()]))
end

function print_bus_3ph_fault_line_trip_schedule()
    println("\n" * "="^118)
    println("BUS-FAULT + LINE-TRIP DISTURBANCE SCHEDULE")
    println("="^118)
    @printf("Fault-on   : t = %.3f s, 3φ-G fault at bus %d, Zf = %.6g %+.6gim pu\n",
            FAULT_APPLY_TIME_S, BUS_FAULT_BUS, FAULT_R_PU, FAULT_X_PU)
    @printf("Fault-clear: t = %.3f s, remove bus fault and trip/open line %d--%d\n",
            BUS_FAULT_CLEAR_TIME_S, BUS_FAULT_TRIP_LINE_FROM, BUS_FAULT_TRIP_LINE_TO)
    @printf("Fault duration = %.3f s = %.1f ms\n",
            BUS_FAULT_CLEAR_TIME_S - FAULT_APPLY_TIME_S,
            1000.0 * (BUS_FAULT_CLEAR_TIME_S - FAULT_APPLY_TIME_S))
    println("="^118 * "\n")
end

function print_fault_disturbance_schedule()
    if FAULT_SCENARIO == :line_fault_line_trip
        print_bus_fault_line_trip_schedule()
    elseif FAULT_SCENARIO == :bus_fault_line_trip
        print_bus_3ph_fault_line_trip_schedule()
    else
        error("Unsupported FAULT_SCENARIO=$(FAULT_SCENARIO). Use :line_fault_line_trip or :bus_fault_line_trip.")
    end
end

bus_index_map_for_network_switch(; buses::Vector{Int} = collect(1:39)) =
    Dict{Int,Int}(b => i for (i, b) in enumerate(buses))

function remove_line_from_ybus!(Y::AbstractMatrix{ComplexF64}, fbus::Int, tbus::Int;
                                buses::Vector{Int} = collect(1:39), r::Float64, x::Float64, b::Float64 = 0.0)
    bus_to_idx = bus_index_map_for_network_switch(; buses = buses)
    if !(haskey(bus_to_idx, fbus) && haskey(bus_to_idx, tbus))
        error("Cannot remove line $(fbus)--$(tbus): one or both buses are missing from Ybus bus ordering.")
    end
    z = complex(r, x)
    if abs(z) < 1e-12
        error("Cannot remove line $(fbus)--$(tbus): near-zero impedance r+jx = $(z).")
    end
    y = inv(z)
    ysh = 1im * b / 2.0
    i = bus_to_idx[fbus]
    j = bus_to_idx[tbus]
    Y[i, i] -= y + ysh
    Y[j, j] -= y + ysh
    Y[i, j] += y
    Y[j, i] += y
    return Y
end

function stamp_branch_ybus!(Y::AbstractMatrix{ComplexF64}, fbus::Int, tbus::Int,
                            r::Float64, x::Float64, b::Float64, tap_raw::Float64; buses::Vector{Int})
    bus_to_idx = Dict{Int,Int}(b => i for (i, b) in enumerate(buses))
    if !(haskey(bus_to_idx, fbus) && haskey(bus_to_idx, tbus))
        return Y
    end
    z = complex(r, x)
    if abs(z) < 1e-12
        @warn "Skipping near-zero impedance branch while building Ybus" fbus=fbus tbus=tbus r=r x=x
        return Y
    end
    tap = abs(tap_raw) < 1e-12 ? 1.0 : tap_raw
    y = inv(z)
    ysh = 1im * b / 2.0
    i = bus_to_idx[fbus]
    j = bus_to_idx[tbus]
    Y[i, i] += (y + ysh) / (tap * tap)
    Y[i, j] -= y / tap
    Y[j, i] -= y / tap
    Y[j, j] += y + ysh
    return Y
end

function build_ybus_dense_from_ieee39_fallback(; buses::Vector{Int} = NETWORK_PQ_ALL_BUSES)
    Y = zeros(ComplexF64, length(buses), length(buses))
    for (f, t, r, x, b, tap) in VPD_BRANCH_TABLE
        stamp_branch_ybus!(Y, f, t, r, x, b, tap; buses = buses)
    end
    return Y
end

_branch_key(fbus::Int, tbus::Int) = fbus < tbus ? (fbus, tbus) : (tbus, fbus)

function find_line_parameters(sys::System, fbus::Int, tbus::Int)
    target = _branch_key(fbus, tbus)
    for line in get_components(Line, sys)
        try
            f, t = line_pair(line)
            if _branch_key(f, t) == target
                r = _line_parameter_or_default(line, :get_r, NaN)
                x = _line_parameter_or_default(line, :get_x, NaN)
                b = _line_parameter_or_default(line, :get_b, 0.0)
                if !isfinite(r) || !isfinite(x)
                    error("Found line $(fbus)--$(tbus), but could not read finite r/x values.")
                end
                return (r = r, x = x, b = b, name = PSY.get_name(line))
            end
        catch
        end
    end
    if haskey(IEEE39_BRANCH_PARAM_FALLBACK, target)
        r, x, b = IEEE39_BRANCH_PARAM_FALLBACK[target]
        @warn "Line component not found in sys; using IEEE-39 MATPOWER fallback branch parameters" fbus=fbus tbus=tbus r=r x=x b=b
        return (r = r, x = x, b = b, name = "MATPOWER_branch_$(fbus)_$(tbus)")
    end
    error("Could not find branch parameters for branch $(fbus)--$(tbus).")
end

function stamp_series_branch_by_indices!(Y::AbstractMatrix{ComplexF64}, i::Int, j::Int;
                                         r::Float64, x::Float64, b::Float64 = 0.0)
    z = complex(r, x)
    if abs(z) < 1e-12
        error("Cannot stamp branch: near-zero impedance r+jx = $(z).")
    end
    y = inv(z)
    ysh = 1im * b / 2.0
    Y[i, i] += y + ysh
    Y[j, j] += y + ysh
    Y[i, j] -= y
    Y[j, i] -= y
    return Y
end

function build_line_fault_ybus_kron_reduced(Ybus_prefault::AbstractMatrix{ComplexF64},
                                            line_from::Int, line_to::Int;
                                            fault_location_from::Float64 = FAULT_LOCATION_FROM,
                                            fault_r::Float64 = FAULT_R_PU, fault_x::Float64 = FAULT_X_PU,
                                            buses::Vector{Int} = collect(1:39),
                                            r::Float64, x::Float64, b::Float64 = 0.0)
    if !(0.0 < fault_location_from < 1.0)
        error("FAULT_LOCATION_FROM must be strictly between 0 and 1. Use 0.5 for a midpoint line fault.")
    end
    zf = complex(fault_r, fault_x)
    if abs(zf) < 1e-12
        error("Fault impedance is too close to zero. Use a small nonzero value.")
    end
    yf = inv(zf)
    bus_to_idx = bus_index_map_for_network_switch(; buses = buses)
    if !(haskey(bus_to_idx, line_from) && haskey(bus_to_idx, line_to))
        error("Line-fault buses $(line_from)--$(line_to) are not present in the NetworkSwitch bus ordering.")
    end
    Ybase = copy(Ybus_prefault)
    remove_line_from_ybus!(Ybase, line_from, line_to; buses = buses, r = r, x = x, b = b)
    n = size(Ybase, 1)
    m = n + 1
    Yaug = zeros(ComplexF64, n + 1, n + 1)
    Yaug[1:n, 1:n] .= Ybase
    i = bus_to_idx[line_from]
    j = bus_to_idx[line_to]
    α = fault_location_from
    r1 = α * r; x1 = α * x; b1 = α * b
    r2 = (1.0 - α) * r; x2 = (1.0 - α) * x; b2 = (1.0 - α) * b
    stamp_series_branch_by_indices!(Yaug, i, m; r = r1, x = x1, b = b1)
    stamp_series_branch_by_indices!(Yaug, m, j; r = r2, x = x2, b = b2)
    Yaug[m, m] += yf
    if abs(Yaug[m, m]) < 1e-12
        error("Cannot Kron-reduce line fault point because Ymm is near zero.")
    end
    Yred = Yaug[1:n, 1:n] .- (Yaug[1:n, m] * transpose(Yaug[m, 1:n])) ./ Yaug[m, m]
    return ComplexF64.(Yred)
end

# Ybus for a NetworkSwitch must include the constant-impedance loads, otherwise
# the switched network loses its demand.  The branch-only builder above is kept
# as is for the fault Ybus construction.
function build_switch_ybus_with_loads(sys::System; buses::Vector{Int} = collect(1:39))
    Y = ComplexF64.(build_ybus_dense_from_system(sys; buses = buses))
    bus_to_idx = Dict{Int,Int}(b => i for (i, b) in enumerate(buses))
    n = 0
    for fa in get_components(FixedAdmittance, sys)
        bnum = PSY.get_number(PSY.get_bus(fa))
        haskey(bus_to_idx, bnum) || continue
        Yval = try ComplexF64(PSY.get_Y(fa)) catch; continue end
        Y[bus_to_idx[bnum], bus_to_idx[bnum]] += Yval
        n += 1
    end
    @info "NetworkSwitch Ybus includes constant-impedance loads/shunts" fixed_admittances_stamped = n
    return Y
end

function make_bus_fault_cleared_by_line_trip_perturbations(sys::System;
                                                           line_from::Int = FAULT_LINE_FROM,
                                                           line_to::Int = FAULT_LINE_TO,
                                                           fault_location_from::Float64 = FAULT_LOCATION_FROM,
                                                           t_fault::Float64 = FAULT_APPLY_TIME_S,
                                                           t_clear::Float64 = FAULT_CLEAR_TIME_S,
                                                           fault_r::Float64 = FAULT_R_PU,
                                                           fault_x::Float64 = FAULT_X_PU,
                                                           buses::Vector{Int} = collect(1:39))
    if t_clear <= t_fault
        error("FAULT_CLEAR_TIME_S must be greater than FAULT_APPLY_TIME_S.")
    end
    Ybus_prefault = build_switch_ybus_with_loads(sys; buses = buses)
    lp = find_line_parameters(sys, line_from, line_to)
    Ybus_fault = build_line_fault_ybus_kron_reduced(Ybus_prefault, line_from, line_to;
        fault_location_from = fault_location_from, fault_r = fault_r, fault_x = fault_x,
        buses = buses, r = lp.r, x = lp.x, b = lp.b)
    Ybus_postfault = copy(Ybus_prefault)
    remove_line_from_ybus!(Ybus_postfault, line_from, line_to; buses = buses, r = lp.r, x = lp.x, b = lp.b)
    zf = complex(fault_r, fault_x)
    yf = inv(zf)
    @info "Prepared line fault cleared by line trip NetworkSwitch perturbations" line_fault = "$(line_from)--$(line_to)" fault_location_from = fault_location_from t_fault = t_fault t_clear = t_clear fault_impedance_pu = zf fault_admittance_pu = yf tripped_line_name = lp.name line_r = lp.r line_x = lp.x line_b = lp.b
    perts = PowerSimulationsDynamics.Perturbation[]
    push!(perts, NetworkSwitch(t_fault, sparse(Ybus_fault)))
    push!(perts, NetworkSwitch(t_clear, sparse(Ybus_postfault)))
    return perts
end


function build_bus_fault_ybus(Ybus_prefault::AbstractMatrix{ComplexF64},
                              fault_bus::Int;
                              fault_r::Float64 = FAULT_R_PU,
                              fault_x::Float64 = FAULT_X_PU,
                              buses::Vector{Int} = collect(1:39))
    zf = complex(fault_r, fault_x)
    if abs(zf) < 1e-12
        error("Fault impedance is too close to zero. Use a small nonzero value.")
    end
    bus_to_idx = bus_index_map_for_network_switch(; buses = buses)
    if !haskey(bus_to_idx, fault_bus)
        error("Bus fault target bus $(fault_bus) is not present in the NetworkSwitch bus ordering.")
    end
    Yfault = copy(Ybus_prefault)
    idx = bus_to_idx[fault_bus]
    Yfault[idx, idx] += inv(zf)
    return ComplexF64.(Yfault)
end

function make_bus_3ph_fault_cleared_by_line_trip_perturbations(sys::System;
                                                               fault_bus::Int = BUS_FAULT_BUS,
                                                               trip_line_from::Int = BUS_FAULT_TRIP_LINE_FROM,
                                                               trip_line_to::Int = BUS_FAULT_TRIP_LINE_TO,
                                                               t_fault::Float64 = FAULT_APPLY_TIME_S,
                                                               t_clear::Float64 = BUS_FAULT_CLEAR_TIME_S,
                                                               fault_r::Float64 = FAULT_R_PU,
                                                               fault_x::Float64 = FAULT_X_PU,
                                                               buses::Vector{Int} = collect(1:39))
    if t_clear <= t_fault
        error("BUS_FAULT_CLEAR_TIME_S must be greater than FAULT_APPLY_TIME_S.")
    end
    Ybus_prefault = build_switch_ybus_with_loads(sys; buses = buses)
    lp = find_line_parameters(sys, trip_line_from, trip_line_to)

    Ybus_fault = build_bus_fault_ybus(Ybus_prefault, fault_bus;
        fault_r = fault_r, fault_x = fault_x, buses = buses)

    Ybus_postfault = copy(Ybus_prefault)
    remove_line_from_ybus!(Ybus_postfault, trip_line_from, trip_line_to;
        buses = buses, r = lp.r, x = lp.x, b = lp.b)

    zf = complex(fault_r, fault_x)
    yf = inv(zf)
    @info "Prepared bus 3φ-G fault cleared by line trip NetworkSwitch perturbations" fault_bus = fault_bus tripped_line = "$(trip_line_from)--$(trip_line_to)" t_fault = t_fault t_clear = t_clear fault_duration_s = (t_clear - t_fault) fault_impedance_pu = zf fault_admittance_pu = yf tripped_line_name = lp.name line_r = lp.r line_x = lp.x line_b = lp.b

    perts = PowerSimulationsDynamics.Perturbation[]
    push!(perts, NetworkSwitch(t_fault, sparse(Ybus_fault)))
    push!(perts, NetworkSwitch(t_clear, sparse(Ybus_postfault)))
    return perts
end

function make_fault_perturbations(sys::System)
    if FAULT_SCENARIO == :line_fault_line_trip
        return make_bus_fault_cleared_by_line_trip_perturbations(sys)
    elseif FAULT_SCENARIO == :bus_fault_line_trip
        return make_bus_3ph_fault_cleared_by_line_trip_perturbations(sys)
    else
        error("Unsupported FAULT_SCENARIO=$(FAULT_SCENARIO). Use :line_fault_line_trip or :bus_fault_line_trip.")
    end
end


# Read the first available state series for a device.
function get_first_state_series(results, dev_name::String, syms::Vector{Symbol})
    for s in syms
        try
            t, x = get_state_series(results, (dev_name, s))
            return (t, x, s)
        catch
        end
    end
    return nothing
end


# Rotor-angle export
function unwrap_angle_rad(a::AbstractVector)
    n = length(a)
    n == 0 && return a
    out = copy(a)
    twoπ = 2π
    for k in 2:n
        Δ = out[k] - out[k - 1]
        if Δ > π
            out[k:end] .-= twoπ
        elseif Δ < -π
            out[k:end] .+= twoπ
        end
    end
    return out
end

function first_finite_index(v::AbstractVector)
    for i in eachindex(v)
        isfinite(v[i]) && return i
    end
    return nothing
end

function zero_start(ang::AbstractVector)
    i0 = first_finite_index(ang)
    i0 === nothing && return ang
    return ang .- ang[i0]
end

function uniform_grid_with_event(Tfinal::Float64, dt::Float64, t_event::Float64)
    tg = collect(0.0:dt:Tfinal)
    if !(any(abs.(tg .- t_event) .< 1e-12))
        push!(tg, t_event)
        sort!(tg)
    end
    return tg
end

function interp1_linear(t::AbstractVector, y::AbstractVector, tq::AbstractVector)
    n = length(t)
    n == 0 && return fill(NaN, length(tq))
    out = fill(NaN, length(tq))
    ord = sortperm(t)
    tt = t[ord]
    yy = y[ord]
    tmin = tt[1]
    tmax = tt[end]
    for (k, x) in pairs(tq)
        if x < tmin || x > tmax
            out[k] = NaN
            continue
        end
        i = searchsortedlast(tt, x)
        if i >= length(tt)
            out[k] = yy[end]
        elseif tt[i] == x
            out[k] = yy[i]
        else
            t0 = tt[i]; t1 = tt[i + 1]; y0 = yy[i]; y1 = yy[i + 1]
            α = (x - t0) / (t1 - t0)
            out[k] = (1 - α) * y0 + α * y1
        end
    end
    return out
end

function dyn_inverter_at_bus(sys::System, busnum::Int)
    g = gen_at_bus(sys, busnum)
    target_name = PSY.get_name(g)
    for d in get_components(DynamicInverter, sys)
        if PSY.get_name(d) == target_name
            return d
        end
    end
    error("No DynamicInverter found at bus $busnum (expected name = $target_name)")
end

function first_finite_value(v::AbstractVector)
    for x in v
        isfinite(x) && return x
    end
    return 0.0
end

function grid_with_events(events::Vector{Float64})
    tg = collect(0.0:EXPORT_DT:TF_SIM)
    for tev in events
        if !(any(abs.(tg .- tev) .< 1e-12))
            push!(tg, tev)
        end
    end
    sort!(tg)
    return tg
end

# Network P/Q reconstruction
normalize_pair(i::Int, j::Int) = i < j ? (i, j) : (j, i)

function line_pair(line)
    f = nothing
    t = nothing
    try
        arc = PSY.get_arc(line)
        f = PSY.get_number(arc.from)
        t = PSY.get_number(arc.to)
    catch
        try
            buses = PSY.get_buses(line)
            f = PSY.get_number(buses[1])
            t = PSY.get_number(buses[2])
        catch err
            error("Could not recover terminal buses for line $(PSY.get_name(line)). Error = $err")
        end
    end
    return normalize_pair(Int(f), Int(t))
end

function _complex_vector_or_nothing(y)
    try
        if eltype(y) <: Complex
            return ComplexF64.(y)
        end
    catch
    end
    return nothing
end

function _maybe_extract_complex_voltage_time_value(out)
    if out isa Tuple && length(out) >= 2
        t = Float64.(out[1])
        vc = _complex_vector_or_nothing(out[2])
        vc !== nothing && return t, vc
    end
    if out isa DataFrame
        names_df = names(out)
        t_col = nothing
        for c in ["t_s", "time", "t", "Time"]
            if c in names_df
                t_col = c
                break
            end
        end
        t_col === nothing && return nothing
        for c in ["V", "v", "V_complex", "voltage", "voltage_complex", "bus_voltage", "phasor", "V_phasor"]
            if c in names_df
                vc = _complex_vector_or_nothing(out[:, c])
                vc !== nothing && return Float64.(out[:, t_col]), vc
            end
        end
        for (cr, ci) in [("Vr", "Vi"), ("V_re", "V_im"), ("real", "imag"),
                         ("V_real", "V_imag"), ("voltage_real", "voltage_imag"),
                         ("re", "im"), ("d", "q"), ("vd", "vq"), ("Vd", "Vq")]
            if cr in names_df && ci in names_df
                return Float64.(out[:, t_col]), ComplexF64.(Float64.(out[:, cr]), Float64.(out[:, ci]))
            end
        end
    end
    return nothing
end

function bus_voltage_complex_phasor_series(results, sys::System, bus::Int)
    b = bus_component(sys, bus)
    bus_name = PSY.get_name(b)
    fn_candidates = Symbol[:get_voltage_series, :get_bus_voltage_series, :get_voltage_phasor_series,
        :get_bus_voltage_phasor_series, :get_complex_voltage_series, :get_bus_complex_voltage_series]
    arg_candidates = Any[bus, bus_name, b, (bus,), (bus_name,)]
    for fn_sym in fn_candidates
        fn = _find_callable(fn_sym)
        fn === nothing && continue
        for arg in arg_candidates
            try
                out = fn(results, arg)
                tv = _maybe_extract_complex_voltage_time_value(out)
                if tv !== nothing
                    t, Vc = tv
                    return (t, Vc, @sprintf("Bus %d exact complex bus-voltage phasor from %s", bus, string(fn_sym)), fn_sym, true)
                end
            catch
            end
        end
    end
    try
        vout = voltage_magnitude_pu(results, sys, bus)
        if vout !== nothing
            tV, Vmag, _, vsrc = vout
            angle_fn = _find_callable(:get_voltage_angle_series)
            if angle_fn !== nothing
                for arg in arg_candidates
                    try
                        outθ = angle_fn(results, arg)
                        if outθ isa Tuple && length(outθ) >= 2
                            tθ = Float64.(outθ[1])
                            θraw = Float64.(outθ[2])
                            θ_on_tV = interp1_linear(tθ, unwrap_angle_rad(θraw), Float64.(tV))
                            if any(isfinite.(θ_on_tV))
                                θfill = first_finite_value(θ_on_tV)
                                θ_on_tV = [isfinite(x) ? x : θfill for x in θ_on_tV]
                                Vc = ComplexF64.(Float64.(Vmag) .* exp.(1im .* θ_on_tV))
                                return (Float64.(tV), Vc,
                                        @sprintf("Bus %d exact bus-voltage phasor from |V| and get_voltage_angle_series", bus),
                                        Symbol("Vmag_" * string(vsrc) * "__angle_get_voltage_angle_series"), true)
                            end
                        end
                    catch
                    end
                end
            end
        end
    catch err
        @warn "Could not construct exact bus-voltage phasor from magnitude and angle" bus = bus err = string(err)
    end
    @warn "Exact bus-voltage phasor unavailable; network P/Q may need a generator-device angle fallback" bus = bus model = model_label(bus)
    return nothing
end

function _line_parameter_or_default(line, getter_sym::Symbol, default_value::Float64)
    try
        getter = getfield(PSY, getter_sym)
        return Float64(getter(line))
    catch
        return default_value
    end
end

function build_ybus_dense_from_system(sys::System; buses::Vector{Int} = NETWORK_PQ_ALL_BUSES)
    nb = length(buses)
    bus_to_idx = Dict{Int,Int}(b => i for (i, b) in enumerate(buses))
    Y = zeros(ComplexF64, nb, nb)
    n_stamped = 0
    for line in get_components(Line, sys)
        pair = line_pair(line)
        f, t = pair
        if !(haskey(bus_to_idx, f) && haskey(bus_to_idx, t))
            continue
        end
        r = _line_parameter_or_default(line, :get_r, NaN)
        x = _line_parameter_or_default(line, :get_x, NaN)
        b_sh = _line_parameter_or_default(line, :get_b, 0.0)
        if !isfinite(r) || !isfinite(x)
            @warn "Could not read line R/X for Ybus build" line = PSY.get_name(line) pair = pair
            continue
        end
        z = complex(r, x)
        if abs(z) < 1e-12
            @warn "Skipping near-zero impedance line in Ybus build" line = PSY.get_name(line) pair = pair
            continue
        end
        y = inv(z)
        ysh = 1im * b_sh / 2.0
        i = bus_to_idx[f]
        j = bus_to_idx[t]
        Y[i, i] += y + ysh
        Y[j, j] += y + ysh
        Y[i, j] -= y
        Y[j, i] -= y
        n_stamped += 1
    end
    if n_stamped == 0
        @warn "No Line components were stamped into Ybus; using embedded IEEE-39 MATPOWER branch fallback. This avoids the false all-buses-islanded NetworkSwitch error after DynamicBranch conversion."
        return build_ybus_dense_from_ieee39_fallback(; buses = buses)
    end
    return Y
end

# Branch currents come from PSID's Il_R / Il_I states; they are never replaced
# by an algebraic (V_i - V_j)/(r + jx).

function vpd_dynamic_branch_endpoints(br)
    candidates = Any[br]
    for s in (:branch, :static_branch, :wrapped_branch)
        try
            push!(candidates, getproperty(br, s))
        catch
        end
    end
    for obj in candidates
        try
            arc = PSY.get_arc(obj)
            return Int(PSY.get_number(arc.from)), Int(PSY.get_number(arc.to))
        catch
        end
    end
    error("Cannot establish directed from/to buses for DynamicBranch $(PSY.get_name(br)); cannot assign current sign.")
end

function vpd_dynamic_branch_state(results, br, symbol::Symbol)
    name = string(PSY.get_name(br))
    tried = String[]
    for candidate in (name, "Line " * name, "line " * name)
        push!(tried, candidate)
        try
            t, y = PSID.get_state_series(results, (candidate, symbol))
            return Float64.(t), Float64.(y), candidate
        catch
        end
    end
    error("No $(symbol) dynamic current state for $(name) (tried $(tried)). " *
          "Inspect show_states_initial_value(sim) for branch state names. " *
          "It would be misleading to substitute algebraic Ybus current here.")
end

function vpd_static_transformer_ybus(buses::Vector{Int}, dynamic_pairs::Set{Tuple{Int,Int}})
    Y = zeros(ComplexF64, length(buses), length(buses))
    count = 0
    TRANSFORMERS_DYNAMIC && return Y, count     # every branch is a dynamic RL state
    for (f, t, r, x, b, tap) in VPD_BRANCH_TABLE
        if abs(tap) > 1.0e-12
            pair = normalize_pair(f, t)
            pair in dynamic_pairs && error("Transformer $(pair) counted as both dynamic and static")
            effective_r = if r > 1e-12
                ENABLE_LINE_XR_CAP ? max(r, abs(x) / LINE_XR_CAP) : r
            else
                ENABLE_LINE_XR_CAP ? max(LINE_RMIN, abs(x) / LINE_XR_CAP) : LINE_RMIN
            end
            stamp_branch_ybus!(Y, f, t, effective_r, x, b, tap; buses = buses)
            count += 1
        end
    end
    count > 0 || error("Missing IEEE39 transformer branches")
    return Y, count
end

function vpd_network_pq_from_currents(results, sys::System, buses::Vector{Int}, tgrid,
                                      Vactual::Matrix{ComplexF64})
    idx = Dict{Int,Int}(b => k for (k,b) in enumerate(buses))
    inj = zeros(ComplexF64, length(buses), length(tgrid))
    branches = collect(get_components(DynamicBranch, sys))
    isempty(branches) && error("No DynamicBranch instances available; cannot use current-based network P/Q")
    # A missing converted AC line would invalidate the hybrid network sum.
    remaining_lines = collect(get_components(Line, sys))
    isempty(remaining_lines) || error("Static Line components remain after conversion: $(length(remaining_lines)). " *
                                   "Hybrid current calculation cannot silently omit those elements.")
    pairs = Set{Tuple{Int,Int}}()
    used_names = String[]
    for br in branches
        f, to = vpd_dynamic_branch_endpoints(br)
        haskey(idx,f) && haskey(idx,to) || error("DynamicBranch outside VPD bus set: $f -> $to")
        pair = normalize_pair(f, to)
        pair in pairs && error("Parallel dynamic branches at $pair require distinct branch identifiers; refusing to merge")
        push!(pairs, pair)
        tr, ir, name_r = vpd_dynamic_branch_state(results, br, :Il_R)
        ti, ii, name_i = vpd_dynamic_branch_state(results, br, :Il_I)
        name_r == name_i || error("Current components refer to different branch names")
        current = interp1_linear(tr, ir, tgrid) .+ 1im .* interp1_linear(ti, ii, tgrid)
        all(isfinite, real.(current)) && all(isfinite, imag.(current)) ||
            error("Nonfinite dynamic current at $f -> $to")
        # Direction check at the pre-event equilibrium; this checks relative
        # frame, orientation, per-unit base and nominal RL branch impedance.
        inner = try PSY.get_branch(br) catch; try getproperty(br, :branch) catch; br end end
        r = _line_parameter_or_default(inner, :get_r, NaN)
        x = _line_parameter_or_default(inner, :get_x, NaN)
        if isfinite(r) && isfinite(x) && abs(complex(r,x)) > 1.0e-9
            expected = (Vactual[idx[f],1] - Vactual[idx[to],1]) / complex(r,x)
            measured = current[1]
            if abs(expected) > 1.0e-5 && abs(measured-expected) > max(0.2*abs(expected),0.02)
                error("Initial current mismatch on $(name_r) ($f -> $to): " *
                      "state=$measured expected series-RL=$expected; check state orientation/bases")
            end
        end
        inj[idx[f], :] .+= current
        inj[idx[to], :] .-= current
        push!(used_names, string(name_r))
    end
    # Dynamic branches are RL states; any algebraic transformers are stamped from
    # the case data (none when TRANSFORMERS_DYNAMIC).
    if VPD_REQUIRED_DYNAMIC_LINE_COUNT > 0 && length(branches) != VPD_REQUIRED_DYNAMIC_LINE_COUNT
        error("Expected $(VPD_REQUIRED_DYNAMIC_LINE_COUNT) dynamic Line branches; got $(length(branches)). " *
              "Verify topology and model before accepting network P/Q")
    end
    Ystatic, ntransformers = vpd_static_transformer_ybus(buses, pairs)
    inj .+= Ystatic * Vactual
    S = Vactual .* conj.(inj)
    return (P = real.(S), Q = imag.(S), dynamic_count = length(branches),
            static_count = ntransformers,
            source = TRANSFORMERS_DYNAMIC ? "PSID_dynamic_branch_state_Il_R_Il_I_all_branches" :
                                            "PSID_dynamic_branch_state_Il_R_Il_I_plus_case39_algebraic_transformer_Y")
end

# Power each bus feeds into the series network, on the export grid, from the
# exact terminal phasors and the dynamic branch-current states.
function network_bus_pq_dev(results, sys::System, pf_data::PowerFlowData, steps::Vector{RealizedLoadStep}; buses::Vector{Int} = NETWORK_PQ_ALL_BUSES)
    ev_times = unique(sort([s.t_step for s in steps]))
    t_grid = grid_with_events(ev_times)
    Vphysical = zeros(ComplexF64, length(buses), length(t_grid))
    angle_source = String[]
    angle_is_exact_bus_voltage = Bool[]
    for (i, b) in enumerate(buses)
        vc_out = bus_voltage_complex_phasor_series(results, sys, b)
        vc_out === nothing && error("No exact terminal-voltage phasor at bus $b; network P/Q cannot be exported")
        tv, Vc_raw, _, src, exact_flag = vc_out
        Vmag = interp1_linear(tv, abs.(Vc_raw), t_grid)
        θ_on_grid = interp1_linear(tv, unwrap_angle_rad(angle.(Vc_raw)), t_grid)
        Vphysical[i, :] .= Vmag .* exp.(1im .* θ_on_grid)
        push!(angle_source, string(src))
        push!(angle_is_exact_bus_voltage, exact_flag)
    end
    findfirst(==(VPD_REFERENCE_BUS), buses) === nothing &&
        error("Reference bus $VPD_REFERENCE_BUS is not in the network P/Q bus set")
    all(angle_is_exact_bus_voltage) ||
        error("Network P/Q needs the exact terminal phasor at every bus; approximate phasors at " *
              "$(buses[.!angle_is_exact_bus_voltage])")
    pq = vpd_network_pq_from_currents(results, sys, buses, t_grid, Vphysical)
    # deviations are measured from the EMT's own initial point, not from the power flow
    return (t = t_grid, buses = buses, P = pq.P, Q = pq.Q, Pref = pq.P[:, 1], Qref = pq.Q[:, 1],
            V0 = Vphysical[:, 1],
            pq_method = pq.source, dynamic_branch_count = pq.dynamic_count,
            static_transformer_count = pq.static_count,
            angle_source = angle_source, angle_is_exact_bus_voltage = angle_is_exact_bus_voltage)
end

function save_network_pq_csv(results, sys::System, pf_data::PowerFlowData, bus::Int, outdir::String, steps::Vector{RealizedLoadStep}; net_cache = nothing, fname::Union{Nothing,String} = nothing, run_id::Union{Nothing,String} = nothing)
    mkpath(outdir)
    net = net_cache
    net === nothing && (net = network_bus_pq_dev(results, sys, pf_data, steps; buses = NETWORK_PQ_ALL_BUSES))
    if net === nothing
        @warn "No network bus P/Q series for CSV" bus = bus
        return nothing
    end
    idx = findfirst(==(bus), net.buses)
    if idx === nothing
        @warn "Requested bus missing from network P/Q export" bus = bus buses = net.buses
        return nothing
    end
    Pabs = vec(net.P[idx, :])
    Qabs = vec(net.Q[idx, :])
    ΔP_network = Pabs .- net.Pref[idx]
    ΔQ_network = Qabs .- net.Qref[idx]
    fname === nothing && (fname = @sprintf("network_pq_bus%02d_julia.csv", bus))
    path = joinpath(outdir, fname)
    nrow = length(net.t)
    run_value = run_id === nothing ? "legacy_export" : run_id
    CSV.write(path, DataFrame(
        t_s = net.t,
        delta_p_network_pu = ΔP_network, delta_q_network_pu = ΔQ_network,
        p_network_abs_pu = Pabs, q_network_abs_pu = Qabs,
        asset_type = fill(model_label(bus), nrow),
        run_id = fill(run_value, nrow), bus = fill(bus, nrow),
        signal_kind = fill("positive_bus_to_network_injection", nrow),
        pq_method = fill(net.pq_method, nrow)))
    return path
end

function print_network_pq_reconstruction_summary(net)
    net === nothing && return
    n_exact = count(net.angle_is_exact_bus_voltage)
    println("\nNetwork P/Q source: $(net.pq_method)")
    println("  dynamic RL branches: $(net.dynamic_branch_count), algebraic transformers: $(net.static_transformer_count)")
    println("  buses with exact terminal phasors: $n_exact of $(length(net.buses))")
    for (i, b) in enumerate(net.buses)
        net.angle_is_exact_bus_voltage[i] || println("  bus $b uses a fallback angle: $(net.angle_source[i])")
    end
end

#Voltage and angle extraction
function _as_float_vector(y)
    if eltype(y) <: Complex
        return Float64.(abs.(y))
    end
    return Float64.(y)
end

function _maybe_extract_time_value(out)
    if out isa Tuple && length(out) >= 2
        t = Float64.(out[1])
        y = _as_float_vector(out[2])
        return t, y
    end
    if out isa DataFrame
        names_df = names(out)
        t_col = nothing
        for c in ["t_s", "time", "t", "Time"]
            if c in names_df
                t_col = c
                break
            end
        end
        y_col = nothing
        for c in ["V_pu", "Vm_pu", "Vmag_pu", "voltage_pu", "magnitude", "Vm"]
            if c in names_df
                y_col = c
                break
            end
        end
        if t_col !== nothing && y_col !== nothing
            return Float64.(out[:, t_col]), _as_float_vector(out[:, y_col])
        end
    end
    return nothing
end

function _find_callable(fn_sym::Symbol)
    isdefined(Main, fn_sym) && return getfield(Main, fn_sym)
    isdefined(PSID, fn_sym) && return getfield(PSID, fn_sym)
    return nothing
end

function voltage_magnitude_pu(results, sys::System, bus::Int)
    b = bus_component(sys, bus)
    bus_name = PSY.get_name(b)
    fn_candidates = Symbol[:get_voltage_magnitude_series, :get_bus_voltage_magnitude_series,
        :get_voltage_series, :get_bus_voltage_series]
    arg_candidates = Any[bus, bus_name, b, (bus,), (bus_name,)]
    for fn_sym in fn_candidates
        fn = _find_callable(fn_sym)
        fn === nothing && continue
        for arg in arg_candidates
            try
                out = fn(results, arg)
                tv = _maybe_extract_time_value(out)
                if tv !== nothing
                    t, V = tv
                    return (t, V, @sprintf("Bus %d voltage magnitude from %s", bus, string(fn_sym)), fn_sym)
                end
            catch
            end
        end
    end
    g = nothing
    try
        g = gen_at_bus(sys, bus)
    catch
        @warn "No bus-voltage accessor and no generator/device fallback for voltage magnitude" bus = bus bus_name = bus_name
        return nothing
    end
    dev_name = PSY.get_name(g)
    voltage_state_candidates = Symbol[:V_t, :Vt, :v_t, :vt, :V_mag, :Vmag, :Vm, :v_mag, :vm, :V, :v]
    s = get_first_state_series(results, dev_name, voltage_state_candidates)
    if s !== nothing
        t, Vraw, sym = s
        V = _as_float_vector(Vraw)
        return (t, V, @sprintf("Bus %d voltage magnitude from device state %s", bus, string(sym)), sym)
    end
    vd_candidates = Symbol[:vd, :v_d, :Vd, :V_d, :vd_filter, :vr_filter, :vd_pll, :vr_cnv]
    vq_candidates = Symbol[:vq, :v_q, :Vq, :V_q, :vq_filter, :vi_filter, :vq_pll, :vi_cnv]
    vd_state = get_first_state_series(results, dev_name, vd_candidates)
    vq_state = get_first_state_series(results, dev_name, vq_candidates)
    if vd_state !== nothing && vq_state !== nothing
        t_vd, vd_raw, vd_sym = vd_state
        t_vq, vq_raw, vq_sym = vq_state
        vd = Float64.(vd_raw)
        vq_on_vd_time = interp1_linear(t_vq, Float64.(vq_raw), t_vd)
        V = sqrt.(vd .^ 2 .+ vq_on_vd_time .^ 2)
        return (Float64.(t_vd), V,
                @sprintf("Bus %d voltage magnitude reconstructed from %s/%s", bus, string(vd_sym), string(vq_sym)),
                Symbol("Vmag_from_" * string(vd_sym) * "_" * string(vq_sym)))
    end
    @warn "No EMT bus voltage magnitude series found" bus = bus bus_name = bus_name dev = dev_name
    return nothing
end

function voltage_magnitude_for_bus(results, sys::System, bus::Int;
                                   t_step::Float64 = 1.0, for_csv::Bool = false)
    out = voltage_magnitude_pu(results, sys, bus)
    out === nothing && return nothing
    t, V, lab, src = out
    if for_csv && RESAMPLE_FOR_CSV
        tg = uniform_grid_with_event(TF_SIM, EXPORT_DT, t_step)
        V  = interp1_linear(t, V, tg)
        t  = tg
    end
    i0 = first_finite_index(V)
    V0 = (i0 !== nothing) ? V[i0] : 0.0
    ΔV = V .- V0
    return (t, V, ΔV, V0, lab, src)
end

function bus_voltage_angle_dev_rad_for_bus(results, sys::System, bus::Int;
                                           t_step::Float64 = 1.0, for_csv::Bool = false)
    vc_out = bus_voltage_complex_phasor_series(results, sys, bus)
    if vc_out === nothing
        @warn "Skipping bus-voltage angle export: exact complex bus-voltage phasor unavailable" bus = bus
        return nothing
    end
    t, Vc, lab, src, exact_flag = vc_out
    θ = zero_start(unwrap_angle_rad(angle.(Vc)))
    if for_csv && RESAMPLE_FOR_CSV
        tg = uniform_grid_with_event(TF_SIM, EXPORT_DT, t_step)
        θ  = interp1_linear(t, θ, tg)
        t  = tg
    end
    return (Float64.(t), Float64.(θ), lab * " (bus-voltage angle deviation)", src, exact_flag)
end

function prepare_vpd_export_dir!(vpd_csvdir::String)
    # Never delete the folder (OneDrive/Explorer may hold handles).  Files are
    # run-stamped and the manifest selects the current run.
    mkpath(vpd_csvdir)

    # Fail early with a clear error only when the folder is genuinely not
    # writable. This probe does not touch any existing VPD result.
    probe = joinpath(vpd_csvdir, ".vpd_write_probe_$(getpid()).tmp")
    try
        open(probe, "w") do io
            write(io, "VPD export directory write test\n")
        end
    catch err
        error("VPD export directory is not writable: $vpd_csvdir. " *
              "Close programs holding the folder and check OneDrive permissions. " *
              "Original error: $(sprint(showerror, err))")
    finally
        try
            isfile(probe) && rm(probe; force = true)
        catch
        end
    end

    return vpd_csvdir
end

function write_vpd_voltage_csv(results, sys::System, bus::Int, outdir::String,
                               run_id::String; t_step::Float64 = 1.0)
    out = voltage_magnitude_for_bus(results, sys, bus; t_step = t_step, for_csv = true)
    out === nothing && error("No voltage magnitude series for VPD bus $bus.")
    t, V, ΔV, V0, lab, src = out
    path = joinpath(outdir, @sprintf("vpd_voltage_bus%02d_%s_julia.csv", bus, run_id))
    n = length(t)
    CSV.write(path, DataFrame(
        t_s = Float64.(t), V_pu = Float64.(V), delta_V_pu = Float64.(ΔV),
        V_initial_pu = fill(Float64(V0), n), run_id = fill(run_id, n),
        bus = fill(bus, n), signal_kind = fill("voltage_magnitude", n),
        source = fill(string(src), n), label = fill(lab, n)))
    return path
end

#Terminal-voltage angle only; a rotor, droop or PLL angle is
function write_vpd_angle_csv(results, sys::System, bus::Int, outdir::String,
                             run_id::String; t_step::Float64 = 1.0)
    out = bus_voltage_angle_dev_rad_for_bus(results, sys, bus; t_step = t_step, for_csv = true)
    out === nothing && error("No complex bus-voltage phasor for VPD angle at bus $bus. " *
                             "Export terminal network voltage rather than a rotor/PLL state.")
    t, angle_rad, lab, src, exact_flag = out
    exact_flag || error("VPD angle at bus $bus uses an approximate, nonexact voltage phasor: $src")
    kind = "network_bus_voltage_angle"
    path = joinpath(outdir, @sprintf("vpd_angle_bus%02d_%s_julia.csv", bus, run_id))
    n = length(t)
    CSV.write(path, DataFrame(
        t_s = Float64.(t), angle_rad = Float64.(angle_rad),
        run_id = fill(run_id, n), bus = fill(bus, n),
        signal_kind = fill(kind, n), source = fill(string(src), n),
        label = fill(lab, n)))
    return (path = path, kind = kind, source = string(src), label = lab)
end

#Save  EMT voltage magnitude and terminal-voltage angle at GFL targets.
function write_vpd_julia_ground_truth(results, sys::System, bus::Int, outdir::String,
                                      run_id::String, voltage_path::String,
                                      angle_path::String; t_step::Float64 = 0.1)
    bus in VPD_UNKNOWN_GFL_BUSES || error("Ground truth requested for non-GFL bus $bus")
    vtab = CSV.read(voltage_path, DataFrame)
    atab = CSV.read(angle_path, DataFrame)
    t = Float64.(vtab.t_s)
    nrow(atab) == length(t) || error("Voltage/angle sample count mismatch at bus $bus")
    any(abs.(t .- Float64.(atab.t_s)) .> 1e-11) &&
        error("Voltage/angle timestamps differ at bus $bus")
    outfile = joinpath(outdir, @sprintf("vpd_julia_ground_truth_bus%02d_%s.csv", bus, run_id))
    CSV.write(outfile, DataFrame(
        t_s = t,
        V_pu = Float64.(vtab.V_pu),
        bus_voltage_angle_dev_rad = Float64.(atab.angle_rad),
        run_id = fill(run_id, length(t)),
        bus = fill(bus, length(t)),
        signal_kind = fill("EMT_REFERENCE_ONLY_NEVER_VPD_INPUT", length(t)),
    ))
    return basename(outfile)
end

#Plot the ten selected buses. The CSV export still covers all 39 buses.
function plot_vpd_julia_reference_package(outdir::String, manifest_df::DataFrame)
    selected = sort(unique(FIGURE_BUSES))
    length(selected) == 10 || error("FIGURE_BUSES must contain exactly ten distinct buses")
    rows = Dict(Int(r.bus) => r for r in eachrow(manifest_df))
    all(b -> haskey(rows, b), selected) || error("Plot buses missing from export manifest")

    for (column, filename, ylabel, plot_title) in (
        (:V_pu, "emt_voltage_magnitude_selected_buses.png", "Voltage magnitude [pu]",
         "EMT bus-voltage magnitude"),
        (:angle_rad, "emt_bus_voltage_angle_selected_buses.png", "Angle deviation [rad]",
         "EMT bus-voltage angle in the network reference frame"),
    )
        panels = Any[]
        for bus in selected
            row = rows[bus]
            csvname = column == :V_pu ? String(row.voltage_file) : String(row.angle_file)
            tab = CSV.read(joinpath(outdir, csvname), DataFrame)
            p = plot(Float64.(tab.t_s), Float64.(tab[!, column]);
                     linewidth = 1.6, label = "", color = :black,
                     title = "Bus $bus", xlabel = "Time [s]", ylabel = ylabel,
                     legend = false, grid = true)
            for te in fault_event_times()
                vline!(p, [te]; color = :gray, linestyle = :dash,
                       linewidth = 0.8, label = "")
            end
            push!(panels, p)
        end
        p_all = plot(panels...; layout = (5, 2), size = (1200, 1450),
                     plot_title = plot_title)
        savefig(p_all, joinpath(outdir, filename))
    end
    return nothing
end

function export_vpd_validation_package(results, sys::System, pf_data::PowerFlowData,
                                       vpd_csvdir::String,
                                       realized_steps::Vector{RealizedLoadStep}, net;
                                       t_step::Float64 = 1.0)
    sort(VPD_UNKNOWN_GFL_BUSES) == sort(GFL_BUSES) ||
        error("VPD targets must equal Julia's physical GFL buses")
    isempty(intersect(VPD_UNKNOWN_GFL_BUSES, VPD_KNOWN_PHASOR_BUSES)) ||
        error("VPD target/known partition overlap")
    sort(unique(vcat(VPD_UNKNOWN_GFL_BUSES,VPD_KNOWN_PHASOR_BUSES))) == collect(30:39) ||
        error("VPD target/known generator partition does not cover IEEE39 buses 30..39")
    VPD_ALLOW_ESTIMATION_TARGET_VOLTAGE_AS_INPUT && error("Unsafe target voltage input setting")
    net === nothing && error("VPD CSV export requires the network P/Q cache.")
    if any(!flag for flag in net.angle_is_exact_bus_voltage)
        bad = [net.buses[i] for i in eachindex(net.buses) if !net.angle_is_exact_bus_voltage[i]]
        error("Cannot advertise exact VPD EMT P/Q: nonexact network-voltage phasor at buses $(bad). " *
              "Correct the PSID bus-voltage phasor extraction; do not use rotor/PLL fallback.")
    end
    # All exported samples use the same uniform time grid.
    RESAMPLE_FOR_CSV || error("VPD validation export requires RESAMPLE_FOR_CSV=true for a consistent time grid.")
    (isfinite(EXPORT_DT) && EXPORT_DT > 0.0) || error("EXPORT_DT must be finite and positive.")
    (SOLVER_DTMAX <= EXPORT_DT + 1e-12) ||
        error("EMT solver dtmax must not exceed the high-rate VPD output grid")
    EXPORT_DT <= 5e-4 || error("EMT export grid is too coarse for the EM VPD validation")
    @info "VPD export grid" export_dt_s = EXPORT_DT sample_rate_hz = 1.0 / EXPORT_DT
    prepare_vpd_export_dir!(vpd_csvdir)
    run_id = Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")
    # Export the EMT initial operating point alongside the power-flow solution.
    julia_pf_file = @sprintf("vpd_julia_pf_bus_%s.csv", run_id)
    CSV.write(joinpath(vpd_csvdir, julia_pf_file), operating_point_table(pf_data, net))
    network_model_file = write_vpd_network_model_csv(sys, vpd_csvdir, run_id)
    manifest_rows = NamedTuple[]

    # every bus with network P/Q is exported
    sort(VPD_EXPORT_BUSES) == sort(net.buses) ||
        error("VPD_EXPORT_BUSES must equal the network P/Q bus set $(net.buses)")

    for bus in VPD_EXPORT_BUSES
        voltage_path = write_vpd_voltage_csv(results, sys, bus, vpd_csvdir, run_id;
                                             t_step = t_step)
        angle_info = write_vpd_angle_csv(results, sys, bus, vpd_csvdir, run_id;
                                         t_step = t_step)

        pq_name = @sprintf("vpd_network_pq_bus%02d_%s_julia.csv", bus, run_id)
        pq_path = save_network_pq_csv(results, sys, pf_data, bus, vpd_csvdir,
                                      realized_steps; net_cache = net, fname = pq_name,
                                      run_id = run_id)
        pq_path === nothing && error("No network P/Q CSV for VPD bus $bus.")

        truth_file = bus in VPD_UNKNOWN_GFL_BUSES ?
            write_vpd_julia_ground_truth(results, sys, bus, vpd_csvdir, run_id,
                                         voltage_path, angle_info.path;
                                         t_step = t_step) : ""
        role = bus in VPD_UNKNOWN_GFL_BUSES ? "unknown_truth_only_voltage" :
               bus in VPD_KNOWN_PHASOR_BUSES ? "known_voltage_input" : "passive_pq_input"
        push!(manifest_rows, (
            schema_version = VPD_MANIFEST_SCHEMA,
            vpd_role = role,
            network_frame = VPD_NETWORK_FRAME_NOTE,
            julia_ground_truth_file = truth_file,
            julia_pf_file = julia_pf_file,
            network_model_file = network_model_file,
            run_id = run_id,
            bus = bus,
            reference_bus = VPD_REFERENCE_BUS,
            export_dt_s = EXPORT_DT,
            sample_rate_hz = 1.0 / EXPORT_DT,
            sampling_description = "high_rate_dense_output_grid_not_decimated",
            angle_kind = angle_info.kind,
            voltage_file = basename(voltage_path),
            angle_file = basename(angle_info.path),
            network_pq_file = basename(pq_path),
            pq_method = net.pq_method,
            dynamic_branch_count = net.dynamic_branch_count,
            static_transformer_count = net.static_transformer_count,
        ))
    end

    manifest_df = DataFrame(manifest_rows)

    #Save both a dated manifest and the current-run manifest.
    archive_manifest_path = joinpath(
        vpd_csvdir,
        @sprintf("vpd_export_manifest_%s.csv", run_id),
    )
    CSV.write(archive_manifest_path, manifest_df)

    manifest_path = joinpath(vpd_csvdir, VPD_MANIFEST_FILENAME)
    CSV.write(manifest_path, manifest_df)
    plot_vpd_julia_reference_package(vpd_csvdir, manifest_df)

    @info "VPD package written" folder = vpd_csvdir run_id = run_id buses = length(VPD_EXPORT_BUSES)
    return manifest_path
end

#Small-signal stability check
function report_small_signal(sim)
    ss = small_signal_analysis(sim)
    @info "Small-signal analysis" stable = ss.stable
    evs = try ss.eigenvalues catch; nothing end
    evs === nothing && return ss
    worst = nothing; worstζ = Inf
    println("\n" * "-"^60); println("LEAST-DAMPED OSCILLATORY MODES (ζ < 5%)")
    @printf("%10s %10s %8s %8s\n", "real", "imag", "f[Hz]", "ζ[%]"); println("-"^60)
    for λ in evs
        ω = imag(λ); ω == 0 && continue; abs(λ) < 1e-6 && continue
        ζ = -real(λ) / abs(λ)
        ζ < worstζ && (worstζ = ζ; worst = λ)
        ζ < 0.05 && @printf("%10.4f %10.4f %8.3f %8.2f\n", real(λ), ω, ω / (2π), ζ * 100)
    end
    println("-"^60)
    worst !== nothing && @printf("LEAST-DAMPED: f=%.3f Hz  ζ=%.2f%%  -> %s\n",
        imag(worst) / (2π), worstζ * 100, worstζ > 0 ? "stable" : "UNSTABLE")
    println("-"^60)
    unstable = sort([λ for λ in evs if real(λ) > 1e-3]; by = real, rev = true)
    if !isempty(unstable)
        println("\n" * "="^60); println("PARTICIPATION IN UNSTABLE MODES (driving states/devices)")
        println("="^60)
        P = nothing
        try P = summary_participation_factors(ss) catch e; @warn "participation unavailable" err = string(e) end
        if P !== nothing
            cols = names(P)
            for λr in unique(round.(real.(unstable); digits = 2))
                idx = findfirst(e -> isapprox(real(e), λr; atol = 5e-2) && imag(e) >= 0, evs)
                idx === nothing && continue
                e = evs[idx]
                @printf("\nMODE  real=%.2f  f=%.2f Hz:\n", real(e), imag(e) / (2π))
                idx + 1 <= length(cols) || continue
                sub = sort(DataFrame(state = P[!, 1], p = abs.(P[!, cols[idx + 1]])), :p, rev = true)
                for r in eachrow(first(sub, 8)); @printf("    %-42s %.3f\n", string(r.state), r.p); end
            end
        end
        println("="^60)
    end
    return ss
end

function make_virtual_inertia_block(; Ta::Float64 = 0.1, kd::Float64 = 0.0, kω::Float64 = 0.0, P_ref::Float64 = 1.0)
    try
        return VirtualInertia(Ta = Ta, kd = kd, kω = kω, P_ref = P_ref)
    catch
        return VirtualInertia(Ta = Ta, kd = kd, k_w = kω, P_ref = P_ref)
    end
end


# Attach SG dynamic model
function attach_synchronous_generator!(sys::System, g::Generator;
    H::Float64,
    D::Float64 = D_SG_DEFAULT,
    droop_R::Float64 = R_SG_PU,
    machine_params::Dict{Symbol,Float64},
    avr_params::Dict{Symbol,Float64}
)
    machine = RoundRotorQuadratic(
        R = machine_params[:Rs],
        Td0_p = machine_params[:Tdo_p],
        Td0_pp = machine_params[:Tdo_pp],
        Tq0_p = machine_params[:Tqo_p],
        Tq0_pp = machine_params[:Tqo_pp],
        Xd = machine_params[:Xd],
        Xq = machine_params[:Xq],
        Xd_p = machine_params[:Xd_p],
        Xq_p = machine_params[:Xq_p],
        Xd_pp = machine_params[:Xd_pp],
        Xl = machine_params[:Xl],
        Se = SG_MACHINE_SATURATION,
    )

    shaft = SingleMass(H = H, D = D)

    avr = AVRTypeI(
        Ka = avr_params[:Ka],
        Ke = avr_params[:Ke],
        Kf = avr_params[:Kf],
        Ta = avr_params[:Ta],
        Te = avr_params[:Te],
        Tf = avr_params[:Tf],
        Tr = avr_params[:Tr],
        Va_lim = SG_AVR_VOLTAGE_LIMITS,
        Ae = SG_AVR_AE,
        Be = SG_AVR_BE,
    )

    tg = TGTypeI(
        R = droop_R,
        Ts = TG_TS,
        Tc = TG_EPS,
        T3 = TG_EPS,
        T4 = TG_EPS,
        T5 = TG_EPS,
        valve_position_limits = SG_GOVERNOR_VALVE_LIMITS,
    )

    pss = PSSFixed(V_pss = SG_PSS_OUTPUT)

    dyn = DynamicGenerator(
        name = PSY.get_name(g),
        ω_ref = DEVICE_REFERENCE_SPEED,
        machine = machine,
        shaft = shaft,
        avr = avr,
        prime_mover = tg,
        pss = pss,
    )

    add_component!(sys, dyn, g)
    return dyn
end

# Attach GFM inverter dynamic model
function attach_grid_forming_inverter!(sys::System, g::Generator;
    droop_R::Float64 = R_GFM_PU,
    kq::Float64 = GFM_KQ,
    Ta::Float64 = GFM_TA_DEFAULT,
    kd::Float64 = GFM_KD_DEFAULT,
    kω::Float64 = GFM_KW_DEFAULT,
    kpv::Float64 = GFM_KPV,
    kiv::Float64 = GFM_KIV,
    kpc::Float64 = GFM_KPC,
    kic::Float64 = GFM_KIC,
    kffv::Float64 = GFM_KFFV,
    kffi::Float64 = GFM_KFFI,
    rv::Float64 = GFM_RV,
    lv::Float64 = GFM_LV,
    ωad::Float64 = GFM_OMEGA_AD,
    kad::Float64 = GFM_KAD,
    lf::Float64 = FILTER_LF,
    rf::Float64 = FILTER_RF,
    cf::Float64 = FILTER_CF,
    lg::Float64 = FILTER_LG,
    rg::Float64 = FILTER_RG
)
    p_ref0, q_ref0 = get_generator_pref_qref(g)
    v_ref0 = get_bus_vref(sys, g)

    converter = AverageConverter(rated_voltage = CONVERTER_RATED_VOLTAGE, rated_current = CONVERTER_RATED_CURRENT)
    vi = make_virtual_inertia_block(Ta = Ta, kd = kd, kω = kω, P_ref = p_ref0)
    reactive = ReactivePowerDroop(kq = kq, ωf = GFM_REACTIVE_FILTER_OMEGA, V_ref = v_ref0)
    outer = OuterControl(vi, reactive)

    inner = VoltageModeControl(
        kpv = kpv,
        kiv = kiv,
        kffv = kffv,
        rv = rv,
        lv = lv,
        kpc = kpc,
        kic = kic,
        kffi = kffi,
        ωad = ωad,
        kad = kad,
    )

    dc = FixedDCSource(voltage = INVERTER_DC_VOLTAGE)
    pll = FixedFrequency()
    filt = LCLFilter(lf = lf, rf = rf, cf = cf, lg = lg, rg = rg)

    @info "GFM controller references initialized" bus = PSY.get_number(PSY.get_bus(g)) dev = PSY.get_name(g) P_ref = p_ref0 Q_ref = q_ref0 V_ref = v_ref0 droop_R = droop_R Ta = Ta kd = kd kω = kω kpv = kpv kiv = kiv kpc = kpc kic = kic rf = rf rg = rg

    dyn = DynamicInverter(
        name = PSY.get_name(g),
        ω_ref = DEVICE_REFERENCE_SPEED,
        converter = converter,
        outer_control = outer,
        inner_control = inner,
        dc_source = dc,
        freq_estimator = pll,
        filter = filt,
    )

    add_component!(sys, dyn, g)
    return dyn
end

function attach_grid_following_inverter!(sys::System, g::Generator;
    kp_p::Float64 = GFL_KP_P,
    ki_p::Float64 = GFL_KI_P,
    ωz_p::Float64 = GFL_OMEGA_Z,
    kp_q::Float64 = GFL_KP_Q_DEFAULT,
    ki_q::Float64 = GFL_KI_Q_DEFAULT,
    ωf_q::Float64 = GFL_OMEGA_F_Q_DEFAULT,
    kω_droop::Float64 = min(GFL_KW_DEFAULT, GFL_KW_CAP),
    kpc::Float64 = GFL_KPC,
    kic::Float64 = GFL_KIC_DEFAULT,
    kffv::Float64 = GFL_KFFV,
    pll_ωlp::Float64 = GFL_PLL_OMEGA_LP,
    pll_kp::Float64 = GFL_PLL_KP,
    pll_ki::Float64 = GFL_PLL_KI,
    lf::Float64 = FILTER_LF,
    rf::Float64 = FILTER_RF,
    cf::Float64 = FILTER_CF,
    lg::Float64 = FILTER_LG,
    rg::Float64 = FILTER_RG
)
    converter = AverageConverter(rated_voltage = CONVERTER_RATED_VOLTAGE, rated_current = CONVERTER_RATED_CURRENT)

    p_ref0, q_ref0 = get_generator_pref_qref(g)
    v_ref0 = get_bus_vref(sys, g)

    ap = ActivePowerPI(Kp_p = kp_p, Ki_p = ki_p, ωz = ωz_p, P_ref = p_ref0)
    try
        PSY.get_ext(ap)["Kω"] = kω_droop
    catch
    end

    rp = ReactivePowerPI(Kp_q = kp_q, Ki_q = ki_q, ωf = ωf_q, V_ref = v_ref0, Q_ref = q_ref0)

    outer = OuterControl(ap, rp)

    inner = CurrentModeControl(kpc = kpc, kic = kic, kffv = kffv)

    dc = FixedDCSource(voltage = INVERTER_DC_VOLTAGE)

    pll = KauraPLL(ω_lp = pll_ωlp, kp_pll = pll_kp, ki_pll = pll_ki)

    filt = LCLFilter(lf = lf, rf = rf, cf = cf, lg = lg, rg = rg)

    @info "GFL controller references initialized" bus = PSY.get_number(PSY.get_bus(g)) dev = PSY.get_name(g) P_ref = p_ref0 Q_ref = q_ref0 V_ref = v_ref0 Kω = kω_droop kp_p = kp_p ki_p = ki_p ωz_p = ωz_p kp_q = kp_q ki_q = ki_q ωf_q = ωf_q kpc = kpc kic = kic pll_ωlp = pll_ωlp pll_kp = pll_kp pll_ki = pll_ki rf = rf rg = rg

    dyn = DynamicInverter(
        name = PSY.get_name(g),
        ω_ref = DEVICE_REFERENCE_SPEED,
        converter = converter,
        outer_control = outer,
        inner_control = inner,
        dc_source = dc,
        freq_estimator = pll,
        filter = filt
    )

    add_component!(sys, dyn, g)
    return dyn
end

# =============================================================================
# Model and runtime reporting based on assembled PowerSystems components.
# =============================================================================

function vpd_report_value(obj, key::Symbol)
    obj === nothing && return nothing
    getter = Symbol("get_", key)
    if isdefined(PSY, getter)
        f = getfield(PSY, getter)
        try
            return f(obj)
        catch
        end
    end
    try
        return getproperty(obj, key)
    catch
        return nothing
    end
end

function vpd_report_nstates(obj)
    obj === nothing && return nothing
    v = vpd_report_value(obj, :n_states)
    if v isa Integer
        return Int(v)
    end
    names = vpd_report_value(obj, :states)
    return names isa AbstractVector ? length(names) : nothing
end

function vpd_report_state_names(obj)
    names = vpd_report_value(obj, :states)
    return names isa AbstractVector ? join(string.(names), ", ") : "unavailable"
end

function vpd_report_scalar(value)
    if value isa Number || value isa Bool || value isa AbstractString || value isa Symbol
        return repr(value)
    elseif value isa Tuple || value isa NamedTuple
        return length(value) <= 20 ? sprint(show, value) : nothing
    elseif value isa AbstractVector
        return length(value) <= 20 && all(x -> x isa Number || x isa Bool || x isa Symbol, value) ?
            sprint(show, value) : nothing
    end
    return nothing
end

function vpd_report_component(io::IO, label::AbstractString, obj; nested::Bool = false)
    println(io, "  ", label, ": ", obj === nothing ? "not present" : string(typeof(obj)))
    obj === nothing && return
    n = vpd_report_nstates(obj)
    println(io, "    declared_state_count = ", n === nothing ? "unavailable" : string(n))
    println(io, "    state_symbols = ", vpd_report_state_names(obj))
    props = try propertynames(obj) catch; Symbol[] end
    for prop in props
        prop in (:internal, :ext, :states, :n_states, :name) && continue
        value = try getproperty(obj, prop) catch; nothing end
        printable = vpd_report_scalar(value)
        printable === nothing || println(io, "    ", prop, " = ", printable)
    end
    # Include scalar extension metadata (e.g. GFL active-power frequency droop Kω).
    extension = vpd_report_value(obj, :ext)
    if extension isa AbstractDict
        for key in sort(collect(keys(extension)); by=string)
            value = vpd_report_scalar(extension[key])
            value === nothing || println(io, "    ext[", repr(key), "] = ", value)
        end
    end
    # OuterControl is itself a composite; log its two underlying loop models.
    if nested
        for prop in props
            prop in (:internal, :ext, :states, :n_states, :name) && continue
            v = try getproperty(obj, prop) catch; nothing end
            v === nothing && continue
            vpd_report_scalar(v) === nothing || continue
            p = lowercase(string(prop))
            if occursin("active", p) || occursin("reactive", p) || occursin("control", p)
                vpd_report_component(io, "$(label).$(prop)", v)
            end
        end
    end
end

function vpd_report_dynamic_injection(sys::System, bus::Int)
    g = gen_at_bus(sys, bus)
    di = vpd_report_value(g, :dynamic_injector)
    if di !== nothing
        return di
    end
    T = is_sg_bus(bus) ? DynamicGenerator : DynamicInverter
    for item in get_components(T, sys)
        if PSY.get_name(item) == PSY.get_name(g)
            return item
        end
    end
    return nothing
end

function vpd_report_hardware(io::IO)
    println(io, "Julia version: ", VERSION)
    println(io, "Julia threads: ", Threads.nthreads())
    println(io, "CPU model: ", try Sys.CPU_NAME catch; "unavailable" end)
    println(io, "CPU threads: ", try Sys.CPU_THREADS catch; "unavailable" end)
    println(io, "Operating system: ", Sys.KERNEL, "  architecture: ", Sys.ARCH)
    println(io, "Total RAM (GiB): ", try round(Sys.total_memory() / 2.0^30, digits=2) catch; "unavailable" end)
    println(io, "Active Julia environment: ", try Base.active_project() catch; "unavailable" end)
    for (label, mod) in (("PowerSystems.jl", PSY), ("PowerSimulationsDynamics.jl", PSID),
                         ("Sundials.jl", Sundials), ("PowerFlows.jl", PowerFlows))
        println(io, label, " version: ", try string(Base.pkgversion(mod)) catch; "unavailable" end)
    end
end

function vpd_report_simulation_summary(io::IO, obj, label::String)
    println(io, "\n", label, " (native package summary):")
    try
        show(io, MIME("text/plain"), obj)
    catch
        try
            show(io, obj)
        catch e
            println(io, "  summary unavailable: ", sprint(showerror, e))
        end
    end
    println(io)
end


#inspection only: the system simulation itself is unchanged.
#The report distinguishes the entire DAE vector from resource/branch states.
function vpd_report_dae_dimension(io::IO, sim)
    println(io, "DAE system dimensions (queried from the assembled simulation):")
    found = false
    for nm in (:problem, :prob, :dae_problem)
        problem = try getproperty(sim, nm) catch; nothing end
        problem === nothing && continue
        u0 = try getproperty(problem, :u0) catch; nothing end
        if u0 isa AbstractVector
            println(io, "  Full DAE initial-state vector dimension: ", length(u0))
            found = true
        end
        du0 = try getproperty(problem, :du0) catch; nothing end
        du0 isa AbstractVector && println(io, "  DAE derivative vector dimension: ", length(du0))
        dvs = try getproperty(problem, :differential_vars) catch; nothing end
        if dvs isa AbstractVector{Bool}
            println(io, "  Differential variables: ", count(identity, dvs))
            println(io, "  Algebraic variables: ", length(dvs) - count(identity, dvs))
        end
        found && break
    end
    if !found
        for nm in (:inputs, :simulation_inputs)
            inp = try getproperty(sim, nm) catch; nothing end
            inp === nothing && continue
            for label in (:n_states, :n_differential_states, :n_algebraic_states)
                val = vpd_report_value(inp, label)
                val === nothing || println(io, "  ", label, ": ", val)
            end
        end
        println(io, "  DAE state vector length not exposed by this installed PSID version; consult simulation summary below.")
    end
end

function vpd_write_emt_complexity_report(sys::System, sim, results,
                                         path::AbstractString;
                                         time_setup_s::Real, time_integrate_s::Real,
                                         time_read_s::Real, time_reconstruct_s::Real,
                                         time_export_s::Real, time_total_s::Real,
                                         manifest::AbstractString,
                                         n_load_events::Int, n_perturbations::Int)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "THREE-PHASE-FAULT EMT MODEL COMPLEXITY AND COMPUTATIONAL EFFORT")
        println(io, "="^78)
        println(io, "Created: ", Dates.now())
        println(io, "Method: runtime introspection of the assembled PowerSystems/PSID model")
        println(io, "All timing values are measured wall-clock seconds on this machine.")
        println(io, "Timings include Julia execution/compilation costs incurred in this run.")

        println(io, "\n1. EXPERIMENT AND EXCITATION\n", "-"^78)
        println(io, "Scenario: ", DISTURBANCE_SCENARIO)
        println(io, "Scenario description: ", describe_scenario())
        println(io, "System base: ", SYS_BASE_MVA, " MVA")
        println(io, "Fundamental frequency: ", F_BASE, " Hz")
        println(io, "Simulation interval: 0.0 to ", TF_SIM, " s")
        println(io, "Load-change enabled: ", LOAD_CHANGE_ENABLE)
        println(io, "Fault enabled: ", BUS_FAULT_ENABLE)
        println(io, "Scheduled load events: ", n_load_events)
        println(io, "Scheduled switching events: ", length(fault_event_times()))
        println(io, "Fault type: three-phase-to-ground")
        println(io, "Fault specification: ", FAULT_SCENARIO)
        println(io, "Fault branch (from-to): ", FAULT_LINE_FROM, "--", FAULT_LINE_TO)
        println(io, "Fault location fraction from sending end: ", FAULT_LOCATION_FROM)
        println(io, "Fault application: ", FAULT_APPLY_TIME_S, " s")
        println(io, "Fault clearing: ", FAULT_CLEAR_TIME_S, " s")
        println(io, "Fault duration: ", BUS_FAULT_DURATION_S, " s")
        println(io, "Fault resistance Rf: ", FAULT_R_PU, " pu")
        println(io, "Fault reactance Xf: ", FAULT_X_PU, " pu")
        println(io, "Fault magnitude |Zf|: ", abs(complex(FAULT_R_PU, FAULT_X_PU)), " pu")
        println(io, "Fault-on and post-trip topology are imposed by original NetworkSwitch callbacks.")
        println(io, "The line-fault midpoint is reduced via the original Ybus fault construction.")
        println(io, "Perturbation objects passed to PSID: ", n_perturbations)
        println(io, "Sine update interval: ", SINE_UPDATE_DT, " s (configured but inactive in fault-only scenario)")
        println(io, "Export sampling interval: ", EXPORT_DT, " s; nominal rate = ", 1/EXPORT_DT, " Hz")
        println(io, "Resampling enabled: ", RESAMPLE_FOR_CSV)
        println(io, "SG buses: ", SG_BUSES, " | GFM buses: ", GFM_BUSES, " | GFL buses: ", GFL_BUSES)
        println(io, "Known-forming buses for VPD: ", VPD_KNOWN_PHASOR_BUSES)
        println(io, "Reference bus: ", VPD_REFERENCE_BUS)

        println(io, "\n2. NUMERICAL DAE SOLVER\n", "-"^78)
        println(io, "Model formulation: ResidualModel (DAE)")
        println(io, "Integrator: Sundials.IDA, BDF, linear_solver=:Dense")
        println(io, "Maximum BDF order: ", SOLVER_MAX_ORDER)
        println(io, "Absolute tolerance: ", SOLVER_ABSTOL)
        println(io, "Relative tolerance: ", SOLVER_RELTOL)
        println(io, "Initial time step (dtinit): ", SOLVER_INIT_DT, " s")
        println(io, "Maximum internal time step (dtmax): ", SOLVER_DTMAX, " s")
        println(io, "IMPORTANT: dtmax and export interval are NOT fixed solver steps.")
        println(io, "The solver adapts its accepted steps; export sampling is separate.")

        println(io, "\n3. COMPUTATIONAL EFFORT (WALL CLOCK)\n", "-"^78)
        println(io, "Network setup, power flow, devices, perturbations, simulation build and small-signal check: ", round(time_setup_s,digits=6), " s")
        println(io, "IDA execute! integration: ", round(time_integrate_s,digits=6), " s")
        println(io, "read_results: ", round(time_read_s,digits=6), " s")
        println(io, "Network P/Q reconstruction: ", round(time_reconstruct_s,digits=6), " s")
        println(io, "VPD voltage/angle/PQ CSV export and selected-bus figures: ", round(time_export_s,digits=6), " s")
        println(io, "TOTAL main() wall clock through CSV export: ", round(time_total_s,digits=6), " s")
        println(io, "Total excludes writing this report and copying the Julia source snapshot.")
        println(io, "For a cold-start end-to-end benchmark, time the Julia process externally.")

        println(io, "\n4. EXECUTION HARDWARE / SOFTWARE\n", "-"^78)
        vpd_report_hardware(io)

        println(io, "\n5. POWER-SYSTEM RESOURCE MODELS AND COMPONENT STATE COUNTS\n", "-"^78)
        println(io, "Synchronous generators: RoundRotorQuadratic machine, SingleMass shaft,")
        println(io, "  AVRTypeI excitation, TGTypeI prime mover/governor, PSSFixed.")
        println(io, "Grid-forming inverters: AverageConverter, VirtualInertia active control,")
        println(io, "  ReactivePowerDroop, VoltageModeControl voltage/current loops,")
        println(io, "  FixedDCSource, FixedFrequency, LCLFilter.")
        println(io, "Grid-following inverters: AverageConverter, ActivePowerPI and")
        println(io, "  ReactivePowerPI outer loops, CurrentModeControl inner loop,")
        println(io, "  KauraPLL, FixedDCSource, LCLFilter.")
        println(io, "Converter models are averaged; individual semiconductor switching is not modeled.")
        println(io, "Component state orders are obtained from assembled models, NOT preassigned")
        println(io, "  as a generic sixth-order SG, 15th-order GFM, or 20th-order GFL.")
        println(io, "The machine, shaft, excitation, prime mover, active/reactive control,")
        println(io, "  inner controllers, PLL, and LCL filter are enumerated separately below.")
        println(io, "State symbols/counts are reported only where exposed by this PSY installation.")
        sum_device_states = 0
        counted = 0
        for (family, buses) in (("SYNCHRONOUS GENERATORS", SG_BUSES),
                                ("GRID-FORMING INVERTERS", GFM_BUSES),
                                ("GRID-FOLLOWING INVERTERS", GFL_BUSES))
            println(io, "\n", family, " (", length(buses), " units)")
            for bus in buses
                dyn = try vpd_report_dynamic_injection(sys, bus) catch; nothing end
                println(io, "\nBus ", bus, " | ", dyn === nothing ? "dynamic resource unavailable" : string(typeof(dyn)))
                if dyn === nothing
                    println(io, "  ERROR: model could not be resolved from assembled system")
                    continue
                end
                counted += 1
                n = vpd_report_nstates(dyn)
                n === nothing || (sum_device_states += n)
                vpd_report_component(io, "Resource", dyn)
                submodels = is_sg_bus(bus) ?
                    (:machine, :shaft, :avr, :prime_mover, :pss) :
                    (:converter, :outer_control, :inner_control, :dc_source, :freq_estimator, :filter, :limiter)
                sub_sum = 0
                sub_count_complete = true
                for sub in submodels
                    comp = vpd_report_value(dyn, sub)
                    vpd_report_component(io, string(sub), comp; nested = sub == :outer_control)
                    if comp !== nothing
                        nc = vpd_report_nstates(comp)
                        if nc === nothing
                            sub_count_complete = false
                        else
                            sub_sum += nc
                        end
                    end
                end
                if sub_count_complete && n !== nothing
                    println(io, "  Subcomponent state-count sum: ", sub_sum,
                            " | resource state count: ", n,
                            " | equality: ", sub_sum == n)
                else
                    println(io, "  NOTE: Some subcomponent counts unavailable; do not infer their values.")
                end
            end
        end
        println(io, "\nDynamic resources resolved: ", counted, " of ", length(SG_BUSES)+length(GFM_BUSES)+length(GFL_BUSES))
        println(io, "Sum of AVAILABLE resource state counts: ", sum_device_states)
        println(io, "Missing state counts must not be treated as zero.")

        println(io, "\n6. TRANSMISSION NETWORK AND BRANCH DYNAMICS\n", "-"^78)
        branches = collect(get_components(DynamicBranch, sys))
        println(io, "Buses in network: ", length(collect(get_components(Bus,sys))))
        println(io, "Dynamic RL branches assembled: ", length(branches))
        println(io, "Dynamic transformer equivalents requested: ", TRANSFORMERS_DYNAMIC)
        println(io, "Transformer equivalent charging b: ", TRANSFORMER_DYN_B, " pu")
        println(io, "Original network: 34 lines plus 12 transformers (when all 46 conversions succeed).")
        println(io, "Transformer replacements use unit-ratio lines; original tap ratios are not retained.")
        println(io, "Series-resistance floor: ", LINE_RMIN, " pu | X/R cap enabled: ", ENABLE_LINE_XR_CAP)
        total_branch_states = 0
        for (i, br) in enumerate(branches)
            nb = vpd_report_nstates(br)
            nb === nothing || (total_branch_states += nb)
            ends = try vpd_dynamic_branch_endpoints(br) catch; ("?", "?") end
            raw = vpd_report_value(br, :branch)
            inner = raw === nothing ? br : raw
            R = vpd_report_value(inner, :r)
            X = vpd_report_value(inner, :x)
            B = vpd_report_value(inner, :b)
            println(io, "  branch ", i, ": ", ends,
                    " | name=", try PSY.get_name(br) catch; "unavailable" end,
                    " | states=", nb === nothing ? "unavailable" : string(nb),
                    " [", vpd_report_state_names(br), "]",
                    " | r=", R, " x=", X, " b=", B)
        end
        println(io, "Sum of AVAILABLE dynamic RL branch state counts: ", total_branch_states)
        println(io, "Each dynamic RL branch has direct/quadrature currents when fully modeled.")
        println(io, "Resource + RL branch state-count subtotal: ", sum_device_states+total_branch_states)
        println(io, "This subtotal EXCLUDES bus-voltage and any additional DAE variables.")
        println(io, "This subtotal is not a verified total if some component state metadata are unavailable.")
        println(io, "PSID may treat d/q bus voltages as differential variables on dynamic networks.")
        println(io, "The PSID simulation summary below is the authoritative assembled DAE dimension.")
        println(io, "Passive standard loads are converted to FixedAdmittance where applicable.")
        println(io, "LoadChange perturbations in this fault-only run: ", n_load_events)
        println(io, "NetworkSwitch perturbations submitted: ", n_perturbations)
        println(io, "FixedAdmittance elements add no independent differential states.")

        println(io, "\n7. NATIVE SIMULATION AND RESULTS SUMMARIES\n", "-"^78)
        vpd_report_dae_dimension(io, sim)
        vpd_report_simulation_summary(io, sim, "PSID SIMULATION")
        vpd_report_simulation_summary(io, results, "PSID RESULTS")

        println(io, "\n8. PARAMETERS, FILES AND REPRODUCIBILITY\n", "-"^78)
        println(io, "Main simulation source: ", @__FILE__)
        println(io, "Network case generator: case39_matpower.jl")
        println(io, "Physical model and experiment settings: emt_parameters.jl")
        println(io, "Figure buses: ", FIGURE_BUSES)
        println(io, "Export: voltage, angle and network P/Q only; no frequency CSV files.")
        println(io, "Written MATPOWER network: ", joinpath(pwd(), "modified_ieee39.m"))
        println(io, "VPD export manifest: ", manifest)
        println(io, "VPD export folder: ", dirname(path))
        println(io, "Every unit's assembled controller, filter and machine parameters are listed in Section 5.")
        println(io, "Every dynamic branch's r/x/b parameters are listed in Section 6.")
        println(io, "All three source files are copied into this folder after report generation.")
        println(io, "WARNING: Exported VPD post-processing may have its own filters; this report describes EMT internals only.")
    end
    return path
end

##Run the fault case and write voltage, angle and network P/Q exports.
function main(; show_system::Bool = true)
    vpd_main_start_ns = time_ns()
    validate_partitions!(SG_BUSES, GFM_BUSES, GFL_BUSES)

    outdir     = pwd()
    plotsdir   = joinpath(outdir, "plots")
    csvdir     = joinpath(plotsdir, "csv")
    vpd_csvdir = joinpath(csvdir, VPD_EXPORT_SUBDIR)
    mkpath(vpd_csvdir)

    realized_steps = build_realized_load_steps_from_events()
    load_event_times = event_times_from_realized_steps(realized_steps)
    BUS_FAULT_ENABLE && (load_event_times = unique(sort(vcat(load_event_times, fault_event_times()))))
    LOAD_CHANGE_ENABLE && print_load_change_schedule(realized_steps)
    BUS_FAULT_ENABLE && print_fault_disturbance_schedule()
    t_event = minimum(load_event_times; init = 1.0)

    #network and power flow -------------------------------------------------
    casefile = joinpath(outdir, "modified_ieee39.m")
    write_modified_ieee39_case(casefile)
    @info "Wrote MATPOWER case" file = casefile

    sys = System(casefile)
    set_units_base_system!(sys, "SYSTEM_BASE")
    enforce_load_bases!(sys)
    transformers_to_dynamic_lines!(sys)
    enforce_min_line_r!(sys)
    add_dynamic_lines!(sys)

    pf_data = PowerFlowData(ACPowerFlow(), sys)
    solve_powerflow!(pf_data)
    update_system_voltages!(sys, pf_data)
    loads_to_constant_impedance!(sys)

    # Operating-point parameters are defined in emt_parameters.jl.
    gen_mbase_map = GENERATOR_BASE_MVA
    mach = SG_MACHINE_PARAMETERS
    stable_avr = SG_AVR_PARAMETERS

    sync_static_gens_to_pf!(sys, pf_data; buses = collect(30:39))
    print_inverter_reference_summary(sys; gfm_buses = GFM_BUSES, gfl_buses = GFL_BUSES)

    for b in 30:39
        g = gen_at_bus(sys, b)
        try
            PSY.set_base_power!(g, gen_mbase_map[b])
        catch e
            @warn "Could not set generator base power" bus = b err = string(e)
        end
    end
    print_dynamic_reference_summary(sys; buses = collect(30:39))

    # --- dynamic models ---------------------------------------------------------
    for b in SG_BUSES
        g = gen_at_bus(sys, b)
        attach_synchronous_generator!(sys, g;
            H = getH(b), D = getDsg(b), droop_R = R_SG_PU,
            machine_params = mach[b], avr_params = stable_avr)
    end
    for b in GFM_BUSES
        g = gen_at_bus(sys, b)
        attach_grid_forming_inverter!(sys, g;
            droop_R = R_GFM_PU, Ta = getGfmTa(b), kd = getGfmKd(b), kω = getGfmKw(b),
            kpv = GFM_KPV, kiv = GFM_KIV, kpc = GFM_KPC, kic = GFM_KIC, ωad = GFM_OMEGA_AD)
    end
    for b in GFL_BUSES
        g = gen_at_bus(sys, b)
        attach_grid_following_inverter!(sys, g;
            kω_droop = min(getGflKw(b), GFL_KW_CAP),
            kp_q = GFL_KP_Q, ki_q = GFL_KI_Q, ωf_q = GFL_OMEGA_F_Q, kic = GFL_KIC)
    end

    if show_system
        @info "System after attaching dynamic models"
        try
            show(sys)
        catch
        end
    end

    # --- disturbances and simulation --------------------------------------------
    perts = PowerSimulationsDynamics.Perturbation[]
    if LOAD_CHANGE_ENABLE
        ensure_perturbable_loads_all_buses!(sys)
        for b in unique([s.bus for s in realized_steps])
            print_loads_on_bus(sys, b)
        end
        append!(perts, make_load_changes_from_event_schedule(sys, realized_steps))
    end
    BUS_FAULT_ENABLE && append!(perts, make_fault_perturbations(sys))

    sim = Simulation!(ResidualModel, sys, pwd(), (0.0, TF_SIM), perts)
    try
        report_small_signal(sim)
    catch e
        @warn "Small-signal analysis failed; continuing" err = string(e)
    end

    vpd_setup_time_s = (time_ns() - vpd_main_start_ns) / 1e9
    vpd_emt_integration_s = @elapsed begin

        execute!(
            sim,
            IDA(
                linear_solver = :Dense,
                max_order = SOLVER_MAX_ORDER
            );
            abstol = SOLVER_ABSTOL,
            reltol = SOLVER_RELTOL,
            dtinit = SOLVER_INIT_DT,
            dtmax = SOLVER_DTMAX,
            maxiters = 9_000_000
        )

    end
    vpd_read_start_ns = time_ns()
    results = read_results(sim)
    vpd_emt_read_s = (time_ns() - vpd_read_start_ns) / 1e9

    # --- VPD export -------------------------------------------------------------
    vpd_reconstruct_start_ns = time_ns()
    net = network_bus_pq_dev(results, sys, pf_data, realized_steps; buses = NETWORK_PQ_ALL_BUSES)
    net === nothing && error("Network P/Q could not be reconstructed; nothing to export for the VPD")
    print_network_pq_reconstruction_summary(net)
    vpd_reconstruct_s = (time_ns() - vpd_reconstruct_start_ns) / 1e9

    vpd_export_start_ns = time_ns()
    manifest = export_vpd_validation_package(results, sys, pf_data, vpd_csvdir,
                                             realized_steps, net; t_step = t_event)
    vpd_export_s = (time_ns() - vpd_export_start_ns) / 1e9
    vpd_total_main_s = (time_ns() - vpd_main_start_ns) / 1e9

    # Reporting is isolated from the plant simulation and VPD calculations.
    try
        report_path = joinpath(vpd_csvdir, "EMT_Model_Complexity_and_Runtime_Report.txt")
        vpd_write_emt_complexity_report(sys, sim, results, report_path;
            time_setup_s=vpd_setup_time_s,
            time_integrate_s=vpd_emt_integration_s,
            time_read_s=vpd_emt_read_s,
            time_reconstruct_s=vpd_reconstruct_s,
            time_export_s=vpd_export_s,
            time_total_s=vpd_total_main_s,
            manifest=manifest,
            n_load_events=length(realized_steps),
            n_perturbations=length(perts))
        for source_name in ("run_ieee39_fault.jl", "case39_matpower.jl", "emt_parameters.jl")
            cp(joinpath(@__DIR__, source_name), joinpath(vpd_csvdir, source_name); force=true)
        end
        @info "EMT complexity report saved" file=report_path
        @printf("EMT integration wall time: %.3f s | total through CSV export: %.3f s\n",
                vpd_emt_integration_s, vpd_total_main_s)
    catch e
        @warn "EMT reporting failed, but original simulation and CSV export are unaffected" error=sprint(showerror, e)
    end

    @info "Done" vpd_folder = vpd_csvdir manifest = manifest
    return sys, results
end

main(show_system = true)