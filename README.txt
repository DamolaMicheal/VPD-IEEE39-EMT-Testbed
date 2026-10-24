IEEE 39-bus EMT fault study for VPD


The project contains three Julia source files:

  run_ieee39_fault.jl    Model assembly, IDA simulation, voltage/angle/PQ export,
                         two selected-bus plots, runtime and complexity report.
  emt_parameters.jl     Simulation settings, resource/controller/filter parameters,
                         disturbance settings, branch-table verification data.
  case39_matpower.jl    Original modified IEEE 39-bus MATPOWER case writer.

To run: open a Julia environment with PowerSystems.jl, PowerFlows.jl,
PowerSimulationsDynamics.jl, Sundials.jl, Plots.jl, CSV.jl and DataFrames.jl.
Keep all three source files in the same folder. From that folder, run:

    julia --project=. run_ieee39_fault.jl

code creates modified_ieee39.m and writes results under:
    plots/csv/Runtime_vpd_fault_line1_2_Bus_030/

Original CSV filenames preserved:
    vpd_voltage_busXX_<run_id>_julia.csv
    vpd_angle_busXX_<run_id>_julia.csv
    vpd_network_pq_busXX_<run_id>_julia.csv
    vpd_julia_ground_truth_busXX_<run_id>.csv
    vpd_network_model_<run_id>.csv
    vpd_julia_pf_bus_<run_id>.csv
    vpd_export_manifest.csv

The  CSV  contains only voltage and angle columns

selected figures (ten buses: 1, 2, 3, 15, 16, 21, 27, 30, 31, 39):
    emt_voltage_magnitude_selected_buses.png
    emt_bus_voltage_angle_selected_buses.png





