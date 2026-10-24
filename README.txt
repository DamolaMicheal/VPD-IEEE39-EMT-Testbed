IEEE 39-bus EMT fault study
============================

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

The original code creates modified_ieee39.m and writes results under:
    plots/csv/Runtime_vpd_fault_line1_2_Bus_030/

Original CSV filenames preserved:
    vpd_voltage_busXX_<run_id>_julia.csv
    vpd_angle_busXX_<run_id>_julia.csv
    vpd_network_pq_busXX_<run_id>_julia.csv
    vpd_julia_ground_truth_busXX_<run_id>.csv
    vpd_network_model_<run_id>.csv
    vpd_julia_pf_bus_<run_id>.csv
    vpd_export_manifest.csv

The ground-truth CSV now contains only voltage and angle columns, not
frequency columns. No terminal-frequency or device-frequency CSVs are created.
The manifest schema is therefore vpd_v7_voltage_angle_network_pq.

Requested figures (ten buses: 1, 2, 3, 15, 16, 21, 27, 30, 31, 39):
    emt_voltage_magnitude_selected_buses.png
    emt_bus_voltage_angle_selected_buses.png

The solver and fault settings, network construction, and voltage/angle/PQ
calculations are carried over from the supplied code. Numeric defaults moved
to emt_parameters.jl are unchanged. The script was statically reviewed but
cannot be executed here because Julia/PSID is unavailable.

Important interoperability note:
    The existing vpd_four_models Python script expects the original
    vpd_v6_all_bus_network_pq manifest with frequency_file and
    device_frequency_file fields. It will NOT accept this new export without
    corresponding changes to the Python manifest reader and frequency scoring.
    Historical CSV files are not deleted. Only new exports omit frequency.

This project does not assert that the NetworkSwitch fault model is identical to
a physical interrupted RL branch model. The same Ybus switching implementation
as the supplied source is retained, including its existing limitations.
