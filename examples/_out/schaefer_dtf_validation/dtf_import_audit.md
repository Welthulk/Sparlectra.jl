# Schäfer DTF/FOR001 import audit

Standalone post-bus `A` cards in cases B-E are interpreted as variant branch-echo markers; the following branch-like payload line is preserved as a branch echo and is not added to the electrical model.

## Case A
- file: /workspace/Sparlectra.jl/data/DTF/FOR001.DAT
- parsed branch count: 27
- parsed bus count: 13
- parsed transformer control count: 5
- parsed outage/trailing-section count: 2 outages; 0 trailing records
- slack bus: BSTADTS1
- PV buses: 
- PQ buses: DELTA2S1, BETA2 S1, ASTADTS1, NORD  S1, OST   S1, WEILERS1, BURG  S1, GAMMA S1, SUED  S1, ALPHA S1, BETA1 S1, DELTA1S1
- branch r ohm: [0.104, 8.1], finite=true
- branch x ohm: [2.94, 34.7], finite=true
- branch g S: [0, 6.19e-06], finite=true
- branch b S: [-1.84e-05, 0.00152], finite=true

### Transformer controls
- source line 1: BETA1 S1 -> BETA2 S1 parallel=C nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 2: BETA1 S1 -> BETA2 S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 3: BETA1 S1 -> BETA2 S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 4: DELTA1S1 -> DELTA2S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=0 max_tap_step=13 added_voltage_angle_deg=0.0
- source line 5: DELTA1S1 -> DELTA2S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=0 max_tap_step=13 added_voltage_angle_deg=0.0

### Node start voltages
- DELTA2S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BETA2 S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ASTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- NORD  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- OST   S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BSTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- WEILERS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BURG  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- GAMMA S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- SUED  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ALPHA S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- BETA1 S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- DELTA1S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true

### Post-bus/trailing records
- none

## Case B
- file: /workspace/Sparlectra.jl/data/DTF/FOR001B.DAT
- parsed branch count: 27
- parsed bus count: 13
- parsed transformer control count: 5
- parsed outage/trailing-section count: 0 outages; 54 trailing records
- slack bus: BSTADTS1
- PV buses: NORD  S1
- PQ buses: DELTA2S1, BETA2 S1, ASTADTS1, OST   S1, WEILERS1, BURG  S1, GAMMA S1, SUED  S1, ALPHA S1, BETA1 S1, DELTA1S1
- branch r ohm: [0.104, 8.1], finite=true
- branch x ohm: [2.94, 34.7], finite=true
- branch g S: [0, 6.19e-06], finite=true
- branch b S: [-1.84e-05, 0.00152], finite=true

### Transformer controls
- source line 1: BETA1 S1 -> BETA2 S1 parallel=C nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 2: BETA1 S1 -> BETA2 S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 3: BETA1 S1 -> BETA2 S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 4: DELTA1S1 -> DELTA2S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=0 max_tap_step=13 added_voltage_angle_deg=0.0
- source line 5: DELTA1S1 -> DELTA2S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=0 max_tap_step=13 added_voltage_angle_deg=0.0

### Node start voltages
- DELTA2S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BETA2 S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ASTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- NORD  S1: start_kV=230.5 vn_kV=230.0 vm_pu=1.0021739130434784 finite=true
- OST   S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BSTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- WEILERS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BURG  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- GAMMA S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- SUED  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ALPHA S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- BETA1 S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- DELTA1S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true

### Post-bus/trailing records
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  DELTA1S1       1.0119.712E+000     0.000 1.48E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  BETA1 S1       0.6035.303E+000     0.000 8.05E-004     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA1 S1  DELTA1S1       1.1621.001E+001     0.000 1.52E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA1 S1  DELTA1S1       1.1419.820E+000     0.000 1.50E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1C  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA2 S1  BURG  S1       2.5011.087E+001     0.000 7.45E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  ASTADTS1       8.1003.470E+001     0.000 2.67E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   WEILERS1  ASTADTS1       0.6002.940E+000     0.000 3.49E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  WEILERS1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  NORD  S1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  WEILERS1       5.5002.320E+001     0.000 1.60E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  BURG  S1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BETA2 S1  WEILERS1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1C  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo

## Case C
- file: /workspace/Sparlectra.jl/data/DTF/FOR001C.DAT
- parsed branch count: 27
- parsed bus count: 13
- parsed transformer control count: 5
- parsed outage/trailing-section count: 0 outages; 54 trailing records
- slack bus: BSTADTS1
- PV buses: 
- PQ buses: DELTA2S1, BETA2 S1, ASTADTS1, NORD  S1, OST   S1, WEILERS1, BURG  S1, GAMMA S1, SUED  S1, ALPHA S1, BETA1 S1, DELTA1S1
- branch r ohm: [0.104, 8.1], finite=true
- branch x ohm: [2.94, 34.7], finite=true
- branch g S: [0, 6.19e-06], finite=true
- branch b S: [-1.84e-05, 0.00152], finite=true

### Transformer controls
- source line 1: BETA1 S1 -> BETA2 S1 parallel=C nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 2: BETA1 S1 -> BETA2 S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 3: BETA1 S1 -> BETA2 S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 4: DELTA1S1 -> DELTA2S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=7 max_tap_step=13 added_voltage_angle_deg=0.0
- source line 5: DELTA1S1 -> DELTA2S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=7 max_tap_step=13 added_voltage_angle_deg=0.0

### Node start voltages
- DELTA2S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BETA2 S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ASTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- NORD  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- OST   S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BSTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- WEILERS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BURG  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- GAMMA S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- SUED  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ALPHA S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- BETA1 S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- DELTA1S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true

### Post-bus/trailing records
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  DELTA1S1       1.0119.712E+000     0.000 1.48E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  BETA1 S1       0.6035.303E+000     0.000 8.05E-004     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA1 S1  DELTA1S1       1.1621.001E+001     0.000 1.52E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA1 S1  DELTA1S1       1.1419.820E+000     0.000 1.50E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1C  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA2 S1  BURG  S1       2.5011.087E+001     0.000 7.45E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  ASTADTS1       8.1003.470E+001     0.000 2.67E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   WEILERS1  ASTADTS1       0.6002.940E+000     0.000 3.49E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  WEILERS1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  NORD  S1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  WEILERS1       5.5002.320E+001     0.000 1.60E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  BURG  S1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BETA2 S1  WEILERS1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1C  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo

## Case D
- file: /workspace/Sparlectra.jl/data/DTF/FOR001D.DAT
- parsed branch count: 27
- parsed bus count: 13
- parsed transformer control count: 5
- parsed outage/trailing-section count: 0 outages; 54 trailing records
- slack bus: S1
- PV buses: 
- PQ buses: DELTA2S1, BETA2 S1, ASTADTS1, NORD  S1, OST   S1, BSTADTS1, WEILERS1, BURG  S1, GAMMA S1, SUED  S1, BETA1 S1, DELTA1S1
- branch r ohm: [0.104, 8.1], finite=true
- branch x ohm: [2.94, 34.7], finite=true
- branch g S: [0, 6.19e-06], finite=true
- branch b S: [-1.84e-05, 0.00152], finite=true

### Transformer controls
- source line 1: BETA1 S1 -> BETA2 S1 parallel=C nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 2: BETA1 S1 -> BETA2 S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 3: BETA1 S1 -> BETA2 S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 4: DELTA1S1 -> DELTA2S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=0 max_tap_step=13 added_voltage_angle_deg=0.0
- source line 5: DELTA1S1 -> DELTA2S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=0 max_tap_step=13 added_voltage_angle_deg=0.0

### Node start voltages
- DELTA2S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BETA2 S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ASTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- NORD  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- OST   S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BSTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- WEILERS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BURG  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- GAMMA S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- SUED  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ALPHA S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- BETA1 S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- DELTA1S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true

### Post-bus/trailing records
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  DELTA1S1       1.0119.712E+000     0.000 1.48E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  BETA1 S1       0.6035.303E+000     0.000 8.05E-004     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA1 S1  DELTA1S1       1.1621.001E+001     0.000 1.52E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA1 S1  DELTA1S1       1.1419.820E+000     0.000 1.50E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1C  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA2 S1  BURG  S1       2.5011.087E+001     0.000 7.45E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  ASTADTS1       8.1003.470E+001     0.000 2.67E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   WEILERS1  ASTADTS1       0.6002.940E+000     0.000 3.49E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  WEILERS1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  NORD  S1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  WEILERS1       5.5002.320E+001     0.000 1.60E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  BURG  S1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BETA2 S1  WEILERS1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1C  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo

## Case E
- file: /workspace/Sparlectra.jl/data/DTF/FOR001E.DAT
- parsed branch count: 27
- parsed bus count: 13
- parsed transformer control count: 5
- parsed outage/trailing-section count: 0 outages; 54 trailing records
- slack bus: BSTADTS1
- PV buses: 
- PQ buses: DELTA2S1, BETA2 S1, ASTADTS1, NORD  S1, OST   S1, WEILERS1, BURG  S1, GAMMA S1, SUED  S1, ALPHA S1, BETA1 S1, DELTA1S1
- branch r ohm: [0.104, 8.1], finite=true
- branch x ohm: [2.94, 34.7], finite=true
- branch g S: [0, 6.19e-06], finite=true
- branch b S: [-1.84e-05, 0.00152], finite=true

### Transformer controls
- source line 1: BETA1 S1 -> BETA2 S1 parallel=C nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 2: BETA1 S1 -> BETA2 S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 3: BETA1 S1 -> BETA2 S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=12.5 actual_tap_step=0 max_tap_step=9 added_voltage_angle_deg=0.0
- source line 4: DELTA1S1 -> DELTA2S1 parallel=B nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=7 max_tap_step=13 added_voltage_angle_deg=60.0
- source line 5: DELTA1S1 -> DELTA2S1 parallel=A nominal_unregulated_kV=400.0 nominal_regulated_kV=231.0 longitudinal_range_percent=18.0 actual_tap_step=7 max_tap_step=13 added_voltage_angle_deg=60.0

### Node start voltages
- DELTA2S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BETA2 S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ASTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- NORD  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- OST   S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BSTADTS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- WEILERS1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- BURG  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- GAMMA S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- SUED  S1: start_kV=230.0 vn_kV=230.0 vm_pu=1.0 finite=true
- ALPHA S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- BETA1 S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true
- DELTA1S1: start_kV=400.0 vn_kV=400.0 vm_pu=1.0 finite=true

### Post-bus/trailing records
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  DELTA1S1       1.0119.712E+000     0.000 1.48E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ALPHA S1  BETA1 S1       0.6035.303E+000     0.000 8.05E-004     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA1 S1  DELTA1S1       1.1621.001E+001     0.000 1.52E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA1 S1  DELTA1S1       1.1419.820E+000     0.000 1.50E-003     3.810"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1C  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BETA2 S1  BURG  S1       2.5011.087E+001     0.000 7.45E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BETA2 S1  BURG  S1       2.5981.087E+001     0.000 7.41E-005     0.725"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BURG  S1  GAMMA S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  OST   S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  ASTADTS1       8.1003.470E+001     0.000 2.67E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   WEILERS1  ASTADTS1       0.6002.940E+000     0.000 3.49E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   DELTA2S1  WEILERS1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   ASTADTS1  NORD  S1       1.2806.360E+000     0.000 7.54E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  ASTADTS1  BSTADTS1       1.3406.720E+000     0.000 1.18E-005     0.700"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  WEILERS1       5.5002.320E+001     0.000 1.60E-004     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BSTADTS1  BURG  S1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1   BETA2 S1  WEILERS1       0.7004.650E+000     0.000 5.50E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1B  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="L1A  BSTADTS1  SUED  S1       1.2305.210E+000     0.000 3.80E-005     1.400"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1C  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  BETA1 S1  BETA2 S1       0.1808.216E+000     0.000-1.84E-005     1.650"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1B  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo
- raw="A"; interpreted kind=variant_branch_echo_marker
- raw="T1A  DELTA1S1  DELTA2S1       0.1048.160E+000     0.000-1.01E-005     2.499"; interpreted kind=variant_branch_echo

