#Parameters for the IEEE 39-bus three-phase-fault EMT study.

const SYS_BASE_MVA = 100.0
const F_BASE       = 60.0
const Ω_BASE_RAD_S = 2π * F_BASE

const TF_SIM              = 0.5
const EXPORT_DT           = 0.000025 
const RESAMPLE_FOR_CSV    = true

const GFM_TA_DEFAULT      = 2.0 * 4.2
const TG_EPS              = 0.001
const TG_TS               = 0.5

# Solver settings
const SOLVER_ABSTOL    = 1e-8
const SOLVER_RELTOL    = 1e-8
const SOLVER_MAX_ORDER = 5
const SOLVER_INIT_DT   = 25e-6
const SOLVER_DTMAX     = 25e-6

const TRANSFORMERS_DYNAMIC = true
const TRANSFORMER_DYN_B    = 0.1    #total charging per converted transformer [pu]


#GFL FFR gains
const GFL_KW_DEFAULT = 75.0
const GFL_KW_map_ext = Dict{Int,Float64}(
    30 => GFL_KW_DEFAULT,
    31 => GFL_KW_DEFAULT,
)

#Zero-resistance branches get LINE_RMIN so every edge has r/x > 0.
const ENABLE_LINE_XR_CAP = false
const LINE_XR_CAP = 25.0
const LINE_RMIN = 0.002

const R_SG_PU  = 0.05
const R_GFM_PU = 0.08

const SG_BUSES  =  [ 33, 35, 37, 38, 39]
const GFM_BUSES =  [31,  36]
const GFL_BUSES = [30, 34, 32]

const GFM_TA_map_ext = Dict{Int,Float64}(
    31 =>  2*3.030,
    32 => 2* 3.580,
    35 =>  3.480,
    34=>  2*4.33,

)

const GFM_KW_DEFAULT = 12.5
const GFM_KW_map_ext = Dict{Int,Float64}(
    31 => GFM_KW_DEFAULT,
    32 => GFM_KW_DEFAULT,
    35 => GFM_KW_DEFAULT,
    38 => GFM_KW_DEFAULT,
)

const GFM_KD_DEFAULT = 12.5
const GFM_KD_map_ext = Dict{Int,Float64}(
    31 => GFM_KD_DEFAULT,
    32 => GFM_KD_DEFAULT,
    35 => GFM_KD_DEFAULT,
    38 => GFM_KD_DEFAULT,
)


#Fault and other disturbance settings.
#The selected case is a midpoint three-phase fault on line 1-2.
const DISTURBANCE_SCENARIO = :line_fault_only

const SCENARIO_SPECS = Dict(
    # Values for the optional load-disturbance experiments
    :sinusoidal_only       => (bus_step = 35, t_step = 0.0, p_step = 0.0,  bus_osc = 35, window = (0.1, 0.5), amp = 0.1, f_hz = 5.0),
    :sinusoidal_then_fault => (bus_step = 35, t_step = 0.0, p_step = 0.0,  bus_osc = 35, window = (0.1, 0.5), amp = 0.1, f_hz = 5.0),
    # Fault-only scenario: no step or sinusoidal load disturbance.
    :line_fault_only       => (bus_step = 15, t_step = 0.0, p_step = 0.0,  bus_osc = 15, window = nothing,    amp = 0.0, f_hz = 0.0),
    :step_then_sinusoidal  => (bus_step = 35, t_step = 0.1, p_step = 2.0,  bus_osc = 35, window = (1.0, 1.5), amp = 0.1, f_hz = 5.0),
    :load_step_only        => (bus_step = 35, t_step = 0.1, p_step = 2.0,  bus_osc = 35, window = nothing,    amp = 0.0, f_hz = 0.0),
    :step_at_target        => (bus_step = 30, t_step = 1.0, p_step = 0.5,  bus_osc = 30, window = nothing,    amp = 0.0, f_hz = 0.0),
    :step_near_target      => (bus_step = 2,  t_step = 1.0, p_step = 0.5,  bus_osc = 2,  window = nothing,    amp = 0.0, f_hz = 0.0),
    :fast_oscillation => (
    bus_step = 27,
    t_step = 0.0,
    p_step = 0.0,
    bus_osc = 27,
    window = (0.2, 0.4),
    amp = 0.1,
    f_hz = 5.0),

)
haskey(SCENARIO_SPECS, DISTURBANCE_SCENARIO) || error("Unknown DISTURBANCE_SCENARIO = $DISTURBANCE_SCENARIO")
const SCENARIO = SCENARIO_SPECS[DISTURBANCE_SCENARIO]
const SINE_UPDATE_DT = 0.000025
const LOAD_STEP_Q_PU = 0.0
const SINE_Q_PU      = 0.0

const LOAD_CHANGE_ENABLE = false
const BUS_FAULT_ENABLE   = DISTURBANCE_SCENARIO == :line_fault_only
const FAULT_SCENARIO     = :line_fault_line_trip   #:line_fault_line_trip or :bus_fault_line_trip

#fault timing and location (only used when BUS_FAULT_ENABLE is true)
const FAULT_APPLY_TIME_S    = 0.10
const BUS_FAULT_DURATION_S      = 0.10

const BUS_FAULT_BUS             = 21
const BUS_FAULT_TRIP_LINE_FROM  = 21
const BUS_FAULT_TRIP_LINE_TO    = 22
const BUS_FAULT_CLEAR_TIME_S    = FAULT_APPLY_TIME_S + BUS_FAULT_DURATION_S
const FAULT_R_PU            = 0.0001
const FAULT_X_PU            = 0.0

#Midpoint fault on line 1-2; cleared by opening the faulted line.
const FAULT_LINE_FROM       = 1
const FAULT_LINE_TO         = 2
const FAULT_LOCATION_FROM   = 0.5
const FAULT_CLEAR_TIME_S    = FAULT_APPLY_TIME_S + BUS_FAULT_DURATION_S

const NETWORK_PQ_ALL_BUSES = collect(1:39)

#Voltage, bus-angle and network P/Q export settings.
const VPD_EXPORT_SUBDIR = "Runtime_vpd_fault_line1_2_Bus_030"
const VPD_MANIFEST_FILENAME = "vpd_export_manifest.csv"
const VPD_MANIFEST_SCHEMA = "vpd_v7_voltage_angle_network_pq"
const VPD_REFERENCE_BUS = 39
const VPD_UNKNOWN_GFL_BUSES = sort(copy(GFL_BUSES))
const VPD_KNOWN_PHASOR_BUSES = sort(unique(vcat(SG_BUSES, GFM_BUSES)))
const VPD_ALLOW_ESTIMATION_TARGET_VOLTAGE_AS_INPUT = false
const VPD_EXPORT_BUSES = collect(1:39)
const VPD_NETWORK_FRAME_NOTE = "psid_network_frame_follows_reference_machine; terminal_voltage_angle_deviation"

#Specific and default Inertia and damping
const H_ALL_map_ext = Dict{Int,Float64}(
    30 => 4.200,
    31 => 3.030,
    32 => 3.580,
    33 => 2.860,
    34 => 4.333,
    35 => 3.480,
    36 => 3.480,
    37 => 2.430,
    38 => 3.450,
    39 => 5.000
)

const D_SCALE = 1.0

#Specifice and default SG damping map
const D_SG_DEFAULT = 40.0 / D_SCALE
const D_SG_map_ext = Dict{Int,Float64}(
    33 => 40.0 / D_SCALE,
    36 => 40.0 / D_SCALE,
    37 => 40.0 / D_SCALE,
    38 => 40.0 / D_SCALE,
    39 => 48.0 / D_SCALE
)


#Synchronous-generator parameters used by the assembled models.
const GENERATOR_BASE_MVA = Dict{Int,Float64}(
    b => (b == 39 ? 10000.0 : b == 34 ? 600.0 : 1000.0)
    for b in 30:39
)
const SG_MACHINE_PARAMETERS = Dict{Int,Dict{Symbol,Float64}}(
    b => Dict{Symbol,Float64}(
        :Rs => 0.0025, :Xl => 0.20,
        :Xd => 2.0, :Xq => 1.9, :Xd_p => 0.40, :Xq_p => 0.60,
        :Xd_pp => 0.30, :Tdo_p => 7.0, :Tqo_p => 0.70,
        :Tdo_pp => 0.035, :Tqo_pp => 0.070,
    ) for b in 30:39
)
const SG_AVR_PARAMETERS = Dict{Symbol,Float64}(
    :Ka => 15.0, :Ke => 1.0, :Kf => 0.08, :Ta => 0.03,
    :Te => 0.8, :Tf => 1.0, :Tr => 0.08,
)
const SG_MACHINE_SATURATION = (0.0, 0.0)
const SG_AVR_VOLTAGE_LIMITS = (-10.0, 10.0)
const SG_AVR_AE = 0.001
const SG_AVR_BE = 1.0
const SG_GOVERNOR_VALVE_LIMITS = (min = -5.0, max = 5.0)
const SG_PSS_OUTPUT = 0.0
const DEVICE_REFERENCE_SPEED = 1.0
const CONVERTER_RATED_VOLTAGE = 1.0
const CONVERTER_RATED_CURRENT = 1.0
const INVERTER_DC_VOLTAGE = 1200.0

#GFM controller and filter parameters.
const GFM_KQ = 0.05
const GFM_REACTIVE_FILTER_OMEGA = 10.0
const GFM_KPV = 1.0
const GFM_KIV = 1000.0
const GFM_KPC = 3.0
const GFM_KIC = 600.0
const GFM_KFFV = 0.0
const GFM_KFFI = 0.0
const GFM_RV = 0.0
const GFM_LV = 0.0
const GFM_OMEGA_AD = 2000.0
const GFM_KAD = 10.0

#GFL power control, current control and PLL parameters.
const GFL_KP_P = 0.8
const GFL_KI_P = 8.0
const GFL_OMEGA_Z = 80.0
const GFL_KP_Q_DEFAULT = 0.8
const GFL_KI_Q_DEFAULT = 12.0
const GFL_OMEGA_F_Q_DEFAULT = 120.0
const GFL_KP_Q = 2.0
const GFL_KI_Q = 30.0
const GFL_OMEGA_F_Q = 500.0
const GFL_KPC = 1.2
const GFL_KIC_DEFAULT = 120.0
const GFL_KIC = 300.0
const GFL_KFFV = 0.0
const GFL_PLL_OMEGA_LP = 80.0
const GFL_PLL_KP = 0.03
const GFL_PLL_KI = 1.5
const GFL_KW_CAP = 35.0

#Six-state LCL filters, shared by GFM and GFL units.
const FILTER_LF = 0.08
const FILTER_RF = 0.05
const FILTER_CF = 0.074
const FILTER_LG = 0.2
const FILTER_RG = 0.04

#Buses selected for voltage and angle figures. All 39 buses remain in the CSV export.
const FIGURE_BUSES = [1, 2, 3, 15, 16, 21, 27, 30, 31, 39]

#Branch parameters used to construct and verify the fault-on and post-trip Y-bus matrices.
#The two modified 16-21 and 21-22 resistances are retained.
const IEEE39_BRANCH_PARAM_FALLBACK = Dict{Tuple{Int,Int},Tuple{Float64,Float64,Float64}}(
    (1, 2) => (0.004, 0.041, 0.698), (1, 39) => (0.001, 0.025, 0.750),
    (2, 3) => (0.001, 0.015, 0.258), (2, 25) => (0.007, 0.009, 0.146),
    (3, 4) => (0.001, 0.021, 0.222), (3, 18) => (0.001, 0.013, 0.214),
    (4, 5) => (0.001, 0.013, 0.134), (4, 14) => (0.001, 0.013, 0.138),
    (5, 6) => (0.000, 0.003, 0.044), (5, 8) => (0.001, 0.011, 0.148),
    (6, 7) => (0.001, 0.009, 0.114), (6, 11) => (0.001, 0.008, 0.138),
    (7, 8) => (0.000, 0.005, 0.078), (8, 9) => (0.002, 0.036, 0.380),
    (9, 39) => (0.001, 0.025, 1.200), (10, 11) => (0.000, 0.004, 0.072),
    (10, 13) => (0.000, 0.004, 0.072), (13, 14) => (0.001, 0.010, 0.172),
    (14, 15) => (0.002, 0.022, 0.366), (15, 16) => (0.001, 0.009, 0.172),
    (16, 17) => (0.001, 0.009, 0.134), (16, 19) => (0.002, 0.020, 0.304),
    (16, 21) => (0.001, 0.014, 0.254), (16, 24) => (0.000, 0.006, 0.068),
    (17, 18) => (0.001, 0.008, 0.132), (17, 27) => (0.001, 0.017, 0.322),
    (21, 22) => (0.001, 0.014, 0.256), (22, 23) => (0.001, 0.010, 0.184),
    (23, 24) => (0.002, 0.035, 0.362), (25, 26) => (0.003, 0.032, 0.512),
    (26, 27) => (0.001, 0.015, 0.240), (26, 28) => (0.004, 0.047, 0.780),
    (26, 29) => (0.006, 0.063, 1.028), (28, 29) => (0.001, 0.015, 0.250),
    (2, 30) => (0.000, 0.018, 0.000), (6, 31) => (0.000, 0.025, 0.000),
    (10, 32) => (0.000, 0.020, 0.000), (12, 11) => (0.002, 0.043, 0.000),
    (12, 13) => (0.002, 0.043, 0.000), (19, 20) => (0.001, 0.014, 0.000),
    (19, 33) => (0.001, 0.014, 0.000), (20, 34) => (0.001, 0.018, 0.000),
    (22, 35) => (0.000, 0.014, 0.000), (23, 36) => (0.000, 0.027, 0.000),
    (25, 37) => (0.001, 0.023, 0.000), (29, 38) => (0.001, 0.016, 0.000),
)

const IEEE39_BRANCH_YBUS_FALLBACK = Tuple{Int,Int,Float64,Float64,Float64,Float64}[
    (1, 2, 0.004, 0.041, 0.698, 0.0), (1, 39, 0.001, 0.025, 0.750, 0.0),
    (2, 3, 0.001, 0.015, 0.258, 0.0), (2, 25, 0.007, 0.009, 0.146, 0.0),
    (3, 4, 0.001, 0.021, 0.222, 0.0), (3, 18, 0.001, 0.013, 0.214, 0.0),
    (4, 5, 0.001, 0.013, 0.134, 0.0), (4, 14, 0.001, 0.013, 0.138, 0.0),
    (5, 6, 0.000, 0.003, 0.044, 0.0), (5, 8, 0.001, 0.011, 0.148, 0.0),
    (6, 7, 0.001, 0.009, 0.114, 0.0), (6, 11, 0.001, 0.008, 0.138, 0.0),
    (7, 8, 0.000, 0.005, 0.078, 0.0), (8, 9, 0.002, 0.036, 0.380, 0.0),
    (9, 39, 0.001, 0.025, 1.200, 0.0), (10, 11, 0.000, 0.004, 0.072, 0.0),
    (10, 13, 0.000, 0.004, 0.072, 0.0), (13, 14, 0.001, 0.010, 0.172, 0.0),
    (14, 15, 0.002, 0.022, 0.366, 0.0), (15, 16, 0.001, 0.009, 0.172, 0.0),
    (16, 17, 0.001, 0.009, 0.134, 0.0), (16, 19, 0.002, 0.020, 0.304, 0.0),

    #Original line 16-21
    # (16, 21, 0.001, 0.014, 0.254, 0.0),
    (16, 21, 0.0005, 0.014, 0.254, 0.0),

    (16, 24, 0.000, 0.006, 0.068, 0.0),
    (17, 18, 0.001, 0.008, 0.132, 0.0), (17, 27, 0.001, 0.017, 0.322, 0.0),

    #Original line 21-22
    #(21, 22, 0.001, 0.014, 0.256, 0.0),
    (21, 22, 0.0005, 0.014, 0.256, 0.0),

    (22, 23, 0.001, 0.010, 0.184, 0.0),
    (23, 24, 0.002, 0.035, 0.362, 0.0), (25, 26, 0.003, 0.032, 0.512, 0.0),
    (26, 27, 0.001, 0.015, 0.240, 0.0), (26, 28, 0.004, 0.047, 0.780, 0.0),
    (26, 29, 0.006, 0.063, 1.028, 0.0), (28, 29, 0.001, 0.015, 0.250, 0.0),
    (2, 30, 0.000, 0.018, 0.000, 1.025), (6, 31, 0.000, 0.025, 0.000, 1.070),
    (10, 32, 0.000, 0.020, 0.000, 1.070), (12, 11, 0.002, 0.043, 0.000, 1.006),
    (12, 13, 0.002, 0.043, 0.000, 1.006), (19, 20, 0.001, 0.014, 0.000, 1.060),
    (19, 33, 0.001, 0.014, 0.000, 1.070), (20, 34, 0.001, 0.018, 0.000, 1.009),
    (22, 35, 0.000, 0.014, 0.000, 1.025), (23, 36, 0.000, 0.027, 0.000, 1.000),
    (25, 37, 0.001, 0.023, 0.000, 1.025), (29, 38, 0.001, 0.016, 0.000, 1.025),
]


#Branch table as the simulated system sees it: with TRANSFORMERS_DYNAMIC the transformer
#rows become unit-ratio lines (tap 0) with charging TRANSFORMER_DYN_B.
const VPD_BRANCH_TABLE = TRANSFORMERS_DYNAMIC ?
    Tuple{Int,Int,Float64,Float64,Float64,Float64}[
        (f, t, r, x, abs(tap) > 1e-12 ? TRANSFORMER_DYN_B : b, 0.0)
        for (f, t, r, x, b, tap) in IEEE39_BRANCH_YBUS_FALLBACK] :
    IEEE39_BRANCH_YBUS_FALLBACK
const VPD_EXPECTED_DYNAMIC_BRANCHES = TRANSFORMERS_DYNAMIC ? 46 : 34
const VPD_EXPECTED_STATIC_TRANSFORMERS = TRANSFORMERS_DYNAMIC ? 0 : 12


const VPD_REQUIRED_DYNAMIC_LINE_COUNT = VPD_EXPECTED_DYNAMIC_BRANCHES
