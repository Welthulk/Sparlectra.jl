# Schäfer DTF/FOR002 validation summary

Generated with native `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!` using executed `transformer_ratio_mode` alternatives. DTF winding/nameplate ratios are preserved for diagnostics; neutral-one mode intentionally does not apply them as neutral off-nominal taps.

## Case A: base case

### Mode neutral_one

- converged: `true`, iterations: `5`, final mismatch: `1.1609699159760757e-13`
- slack: `BSTADTS1` Sparlectra P/Q=309.381/-268.376; FOR002 P/Q=309.4/-268.4; delta=-0.019/0.024
- losses Sparlectra P/Q=7.785/-308.376; FOR002 P/Q=9.4/-308.4; delta=-1.615/0.024
- max deviations: |dV|=0.045 kV (0.0 pu), |dVa|=0.043 deg, branch |dP|=0.197 MW, branch |dQ|=0.045 MVar
- transformer active losses: series=0.113 MW, no-load G/shunt=1.595 MW, total=1.708 MW; total active network loss=7.785 MW
- note: 400-kV mean voltage bias is below 1 kV; old ~2 kV low bias is not visible as a uniform offset.

### Top voltage deviations

- BETA1 S1: model 408.545 kV vs FOR002 408.5 kV; dV=0.045 kV; dVa=-0.024 deg
- DELTA2S1: model 234.432 kV vs FOR002 234.4 kV; dV=0.032 kV; dVa=-0.002 deg
- BETA2 S1: model 234.63 kV vs FOR002 234.6 kV; dV=0.03 kV; dVa=0.041 deg
- OST   S1: model 229.53 kV vs FOR002 229.5 kV; dV=0.03 kV; dVa=0.026 deg
- DELTA1S1: model 409.071 kV vs FOR002 409.1 kV; dV=-0.029 kV; dVa=-0.009 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg
- BURG  S1: model 230.426 kV vs FOR002 230.4 kV; dV=0.026 kV; dVa=0.01 deg
- ASTADTS1: model 231.225 kV vs FOR002 231.2 kV; dV=0.025 kV; dVa=0.007 deg
- WEILERS1: model 234.076 kV vs FOR002 234.1 kV; dV=-0.024 kV; dVa=-0.043 deg
- ALPHA S1: model 403.98 kV vs FOR002 404.0 kV; dV=-0.02 kV; dVa=0.037 deg

### Top branch P deviations

- T1C BETA1 S1 -> BETA2 S1: model -101.897 MW vs FOR002 -101.7 MW; dP=-0.197 MW
- T1A BETA1 S1 -> BETA2 S1: model -101.897 MW vs FOR002 -101.7 MW; dP=-0.197 MW
- T1B BETA1 S1 -> BETA2 S1: model -101.897 MW vs FOR002 -101.7 MW; dP=-0.197 MW
- T1B DELTA1S1 -> DELTA2S1: model -48.189 MW vs FOR002 -48.0 MW; dP=-0.189 MW
- T1A DELTA1S1 -> DELTA2S1: model -48.189 MW vs FOR002 -48.0 MW; dP=-0.189 MW
- L1 DELTA2S1 -> WEILERS1: model -93.446 MW vs FOR002 -93.4 MW; dP=-0.046 MW
- L1 ALPHA S1 -> BETA1 S1: model -270.054 MW vs FOR002 -270.1 MW; dP=0.046 MW
- L1 ALPHA S1 -> DELTA1S1: model -129.946 MW vs FOR002 -129.9 MW; dP=-0.046 MW
- L1 WEILERS1 -> ASTADTS1: model 153.556 MW vs FOR002 153.6 MW; dP=-0.044 MW
- L1B BSTADTS1 -> SUED  S1: model -74.856 MW vs FOR002 -74.9 MW; dP=0.044 MW

### Top branch Q deviations

- L1B BURG  S1 -> GAMMA S1: model 50.445 MVar vs FOR002 50.4 MVar; dQ=0.045 MVar
- L1A BURG  S1 -> GAMMA S1: model 50.445 MVar vs FOR002 50.4 MVar; dQ=0.045 MVar
- L1 WEILERS1 -> ASTADTS1: model 196.955 MVar vs FOR002 197.0 MVar; dQ=-0.045 MVar
- L1 ALPHA S1 -> DELTA1S1: model -96.041 MVar vs FOR002 -96.0 MVar; dQ=-0.041 MVar
- L1 ALPHA S1 -> BETA1 S1: model -103.959 MVar vs FOR002 -104.0 MVar; dQ=0.041 MVar
- L1 BETA2 S1 -> WEILERS1: model 26.56 MVar vs FOR002 26.6 MVar; dQ=-0.04 MVar
- L1 DELTA2S1 -> WEILERS1: model 30.537 MVar vs FOR002 30.5 MVar; dQ=0.037 MVar
- T1B DELTA1S1 -> DELTA2S1: model 23.864 MVar vs FOR002 23.9 MVar; dQ=-0.036 MVar
- T1A DELTA1S1 -> DELTA2S1: model 23.864 MVar vs FOR002 23.9 MVar; dQ=-0.036 MVar
- L1 BSTADTS1 -> WEILERS1: model -44.165 MVar vs FOR002 -44.2 MVar; dQ=0.035 MVar

### Mode winding_over_network

- converged: `true`, iterations: `5`, final mismatch: `1.539879335155092e-13`
- slack: `BSTADTS1` Sparlectra P/Q=309.365/-265.641; FOR002 P/Q=309.4/-268.4; delta=-0.035/2.759
- losses Sparlectra P/Q=7.777/-305.641; FOR002 P/Q=9.4/-308.4; delta=-1.623/2.759
- max deviations: |dV|=1.937 kV (0.005 pu), |dVa|=0.043 deg, branch |dP|=0.24 MW, branch |dQ|=1.169 MVar
- transformer active losses: series=0.112 MW, no-load G/shunt=1.588 MW, total=1.7 MW; total active network loss=7.777 MW
- note: 400-kV mean voltage bias remains visible at -1.894 kV.

### Top voltage deviations

- ALPHA S1: model 402.063 kV vs FOR002 404.0 kV; dV=-1.937 kV; dVa=0.027 deg
- DELTA1S1: model 407.187 kV vs FOR002 409.1 kV; dV=-1.913 kV; dVa=-0.008 deg
- BETA1 S1: model 406.67 kV vs FOR002 408.5 kV; dV=-1.83 kV; dVa=-0.022 deg
- WEILERS1: model 234.048 kV vs FOR002 234.1 kV; dV=-0.052 kV; dVa=-0.042 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg
- GAMMA S1: model 228.474 kV vs FOR002 228.5 kV; dV=-0.026 kV; dVa=0.032 deg
- OST   S1: model 229.522 kV vs FOR002 229.5 kV; dV=0.022 kV; dVa=0.026 deg
- NORD  S1: model 232.881 kV vs FOR002 232.9 kV; dV=-0.019 kV; dVa=0.022 deg
- DELTA2S1: model 234.383 kV vs FOR002 234.4 kV; dV=-0.017 kV; dVa=-0.0 deg
- ASTADTS1: model 231.211 kV vs FOR002 231.2 kV; dV=0.011 kV; dVa=0.008 deg

### Top branch P deviations

- T1B DELTA1S1 -> DELTA2S1: model -48.24 MW vs FOR002 -48.0 MW; dP=-0.24 MW
- T1A DELTA1S1 -> DELTA2S1: model -48.24 MW vs FOR002 -48.0 MW; dP=-0.24 MW
- T1C BETA1 S1 -> BETA2 S1: model -101.865 MW vs FOR002 -101.7 MW; dP=-0.165 MW
- T1A BETA1 S1 -> BETA2 S1: model -101.865 MW vs FOR002 -101.7 MW; dP=-0.165 MW
- T1B BETA1 S1 -> BETA2 S1: model -101.865 MW vs FOR002 -101.7 MW; dP=-0.165 MW
- L1 DELTA2S1 -> WEILERS1: model -93.525 MW vs FOR002 -93.4 MW; dP=-0.125 MW
- L1 BETA2 S1 -> WEILERS1: model 1.079 MW vs FOR002 1.0 MW; dP=0.079 MW
- L1B BETA1 S1 -> DELTA1S1: model 16.926 MW vs FOR002 17.0 MW; dP=-0.074 MW
- L1 ALPHA S1 -> BETA1 S1: model -270.029 MW vs FOR002 -270.1 MW; dP=0.071 MW
- L1 ALPHA S1 -> DELTA1S1: model -129.971 MW vs FOR002 -129.9 MW; dP=-0.071 MW

### Top branch Q deviations

- L1 WEILERS1 -> ASTADTS1: model 195.837 MVar vs FOR002 197.0 MVar; dQ=-1.163 MVar
- L1 BSTADTS1 -> BURG  S1: model -51.274 MVar vs FOR002 -52.4 MVar; dQ=1.126 MVar
- L1 DELTA2S1 -> WEILERS1: model 29.783 MVar vs FOR002 30.5 MVar; dQ=-0.717 MVar
- L1 BETA2 S1 -> WEILERS1: model 25.905 MVar vs FOR002 26.6 MVar; dQ=-0.695 MVar
- T1C BETA1 S1 -> BETA2 S1: model 11.06 MVar vs FOR002 11.7 MVar; dQ=-0.64 MVar
- T1A BETA1 S1 -> BETA2 S1: model 11.06 MVar vs FOR002 11.7 MVar; dQ=-0.64 MVar
- T1B BETA1 S1 -> BETA2 S1: model 11.06 MVar vs FOR002 11.7 MVar; dQ=-0.64 MVar
- L1B BETA1 S1 -> DELTA1S1: model -50.867 MVar vs FOR002 -51.4 MVar; dQ=0.533 MVar
- T1B DELTA1S1 -> DELTA2S1: model 23.37 MVar vs FOR002 23.9 MVar; dQ=-0.53 MVar
- T1A DELTA1S1 -> DELTA2S1: model 23.37 MVar vs FOR002 23.9 MVar; dQ=-0.53 MVar

## Case B: NORD PV bus at 230.5 kV

### Mode neutral_one

- converged: `true`, iterations: `5`, final mismatch: `1.7372965437376723e-13`
- slack: `BSTADTS1` Sparlectra P/Q=309.209/-201.473; FOR002 P/Q=309.2/-201.5; delta=0.009/0.027
- losses Sparlectra P/Q=7.618/-307.993; FOR002 P/Q=9.2/-308.0; delta=-1.582/0.007
- max deviations: |dV|=0.047 kV (0.0 pu), |dVa|=0.05 deg, branch |dP|=0.204 MW, branch |dQ|=0.05 MVar
- transformer active losses: series=0.113 MW, no-load G/shunt=1.591 MW, total=1.703 MW; total active network loss=7.618 MW
- note: NORD effective type=PV, V=230.5 kV, Q=-26.52 MVar.

### Top voltage deviations

- BURG  S1: model 230.247 kV vs FOR002 230.2 kV; dV=0.047 kV; dVa=0.017 deg
- BETA1 S1: model 407.941 kV vs FOR002 407.9 kV; dV=0.041 kV; dVa=-0.008 deg
- ALPHA S1: model 403.341 kV vs FOR002 403.3 kV; dV=0.041 kV; dVa=0.05 deg
- ASTADTS1: model 230.638 kV vs FOR002 230.6 kV; dV=0.038 kV; dVa=0.037 deg
- WEILERS1: model 233.637 kV vs FOR002 233.6 kV; dV=0.037 kV; dVa=-0.02 deg
- OST   S1: model 229.235 kV vs FOR002 229.2 kV; dV=0.035 kV; dVa=0.039 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg
- GAMMA S1: model 228.317 kV vs FOR002 228.3 kV; dV=0.017 kV; dVa=0.038 deg
- BETA2 S1: model 234.316 kV vs FOR002 234.3 kV; dV=0.016 kV; dVa=-0.041 deg
- DELTA2S1: model 234.013 kV vs FOR002 234.0 kV; dV=0.013 kV; dVa=0.017 deg

### Top branch P deviations

- T1C BETA1 S1 -> BETA2 S1: model -101.862 MW vs FOR002 -101.7 MW; dP=-0.162 MW
- T1A BETA1 S1 -> BETA2 S1: model -101.862 MW vs FOR002 -101.7 MW; dP=-0.162 MW
- T1B BETA1 S1 -> BETA2 S1: model -101.862 MW vs FOR002 -101.7 MW; dP=-0.162 MW
- T1B DELTA1S1 -> DELTA2S1: model -48.242 MW vs FOR002 -48.1 MW; dP=-0.142 MW
- T1A DELTA1S1 -> DELTA2S1: model -48.242 MW vs FOR002 -48.1 MW; dP=-0.142 MW
- L1 ASTADTS1 -> NORD  S1: model -99.747 MW vs FOR002 -99.7 MW; dP=-0.047 MW
- L1B BSTADTS1 -> SUED  S1: model -74.856 MW vs FOR002 -74.9 MW; dP=0.044 MW
- L1A BSTADTS1 -> SUED  S1: model -74.856 MW vs FOR002 -74.9 MW; dP=0.044 MW
- L1C BETA2 S1 -> BURG  S1: model 97.756 MW vs FOR002 97.8 MW; dP=-0.044 MW
- L1B BETA2 S1 -> BURG  S1: model 97.756 MW vs FOR002 97.8 MW; dP=-0.044 MW

### Top branch Q deviations

- L1B BURG  S1 -> GAMMA S1: model 50.452 MVar vs FOR002 50.5 MVar; dQ=-0.048 MVar
- L1A BURG  S1 -> GAMMA S1: model 50.452 MVar vs FOR002 50.5 MVar; dQ=-0.048 MVar
- L1 DELTA2S1 -> ASTADTS1: model 16.363 MVar vs FOR002 16.4 MVar; dQ=-0.037 MVar
- L1 ASTADTS1 -> OST   S1: model 46.565 MVar vs FOR002 46.6 MVar; dQ=-0.035 MVar
- L1 ASTADTS1 -> NORD  S1: model 23.735 MVar vs FOR002 23.7 MVar; dQ=0.035 MVar
- L1 WEILERS1 -> ASTADTS1: model 208.468 MVar vs FOR002 208.5 MVar; dQ=-0.032 MVar
- L1C BETA2 S1 -> BURG  S1: model 63.53 MVar vs FOR002 63.5 MVar; dQ=0.03 MVar
- L1B BETA2 S1 -> BURG  S1: model 63.53 MVar vs FOR002 63.5 MVar; dQ=0.03 MVar
- L1B BSTADTS1 -> SUED  S1: model -26.426 MVar vs FOR002 -26.4 MVar; dQ=-0.026 MVar
- L1A BSTADTS1 -> SUED  S1: model -26.426 MVar vs FOR002 -26.4 MVar; dQ=-0.026 MVar

### Mode winding_over_network

- converged: `true`, iterations: `5`, final mismatch: `1.4876988529977098e-13`
- slack: `BSTADTS1` Sparlectra P/Q=309.195/-199.139; FOR002 P/Q=309.2/-201.5; delta=-0.005/2.361
- losses Sparlectra P/Q=7.612/-305.264; FOR002 P/Q=9.2/-308.0; delta=-1.588/2.736
- max deviations: |dV|=1.875 kV (0.005 pu), |dVa|=0.05 deg, branch |dP|=0.236 MW, branch |dQ|=1.218 MVar
- transformer active losses: series=0.113 MW, no-load G/shunt=1.583 MW, total=1.696 MW; total active network loss=7.612 MW
- note: NORD effective type=PV, V=230.5 kV, Q=-26.125 MVar.

### Top voltage deviations

- DELTA1S1: model 406.525 kV vs FOR002 408.4 kV; dV=-1.875 kV; dVa=0.009 deg
- ALPHA S1: model 401.43 kV vs FOR002 403.3 kV; dV=-1.87 kV; dVa=0.039 deg
- BETA1 S1: model 406.072 kV vs FOR002 407.9 kV; dV=-1.828 kV; dVa=-0.006 deg
- DELTA2S1: model 233.967 kV vs FOR002 234.0 kV; dV=-0.033 kV; dVa=0.019 deg
- OST   S1: model 229.229 kV vs FOR002 229.2 kV; dV=0.029 kV; dVa=0.04 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg
- ASTADTS1: model 230.627 kV vs FOR002 230.6 kV; dV=0.027 kV; dVa=0.037 deg
- BURG  S1: model 230.226 kV vs FOR002 230.2 kV; dV=0.026 kV; dVa=0.018 deg
- BETA2 S1: model 234.277 kV vs FOR002 234.3 kV; dV=-0.023 kV; dVa=-0.038 deg
- WEILERS1: model 233.611 kV vs FOR002 233.6 kV; dV=0.011 kV; dVa=-0.019 deg

### Top branch P deviations

- T1B DELTA1S1 -> DELTA2S1: model -48.292 MW vs FOR002 -48.1 MW; dP=-0.192 MW
- T1A DELTA1S1 -> DELTA2S1: model -48.292 MW vs FOR002 -48.1 MW; dP=-0.192 MW
- T1C BETA1 S1 -> BETA2 S1: model -101.83 MW vs FOR002 -101.7 MW; dP=-0.13 MW
- T1A BETA1 S1 -> BETA2 S1: model -101.83 MW vs FOR002 -101.7 MW; dP=-0.13 MW
- T1B BETA1 S1 -> BETA2 S1: model -101.83 MW vs FOR002 -101.7 MW; dP=-0.13 MW
- L1 BETA2 S1 -> WEILERS1: model 0.98 MW vs FOR002 0.9 MW; dP=0.08 MW
- L1 ASTADTS1 -> NORD  S1: model -99.747 MW vs FOR002 -99.7 MW; dP=-0.047 MW
- L1 WEILERS1 -> ASTADTS1: model 153.444 MW vs FOR002 153.4 MW; dP=0.044 MW
- L1B BSTADTS1 -> SUED  S1: model -74.856 MW vs FOR002 -74.9 MW; dP=0.044 MW
- L1A BSTADTS1 -> SUED  S1: model -74.856 MW vs FOR002 -74.9 MW; dP=0.044 MW

### Top branch Q deviations

- L1 WEILERS1 -> ASTADTS1: model 207.282 MVar vs FOR002 208.5 MVar; dQ=-1.218 MVar
- L1 BSTADTS1 -> BURG  S1: model -42.402 MVar vs FOR002 -43.5 MVar; dQ=1.098 MVar
- L1 DELTA2S1 -> WEILERS1: model 30.546 MVar vs FOR002 31.3 MVar; dQ=-0.754 MVar
- L1 BETA2 S1 -> WEILERS1: model 32.187 MVar vs FOR002 32.9 MVar; dQ=-0.713 MVar
- T1C BETA1 S1 -> BETA2 S1: model 10.111 MVar vs FOR002 10.7 MVar; dQ=-0.589 MVar
- T1A BETA1 S1 -> BETA2 S1: model 10.111 MVar vs FOR002 10.7 MVar; dQ=-0.589 MVar
- T1B BETA1 S1 -> BETA2 S1: model 10.111 MVar vs FOR002 10.7 MVar; dQ=-0.589 MVar
- L1A BETA1 S1 -> DELTA1S1: model -49.457 MVar vs FOR002 -50.0 MVar; dQ=0.543 MVar
- L1B BETA1 S1 -> DELTA1S1: model -49.858 MVar vs FOR002 -50.4 MVar; dQ=0.542 MVar
- T1B DELTA1S1 -> DELTA2S1: model 24.323 MVar vs FOR002 24.8 MVar; dQ=-0.477 MVar

## Case C: DELTA fixed tap step 7

### Mode neutral_one

- converged: `true`, iterations: `5`, final mismatch: `7.30526750203353e-14`
- slack: `BSTADTS1` Sparlectra P/Q=311.895/-221.219; FOR002 P/Q=311.6/-223.2; delta=0.295/1.981
- losses Sparlectra P/Q=10.355/-261.219; FOR002 P/Q=11.6/-263.2; delta=-1.245/1.981
- max deviations: |dV|=0.852 kV (0.002 pu), |dVa|=0.049 deg, branch |dP|=1.63 MW, branch |dQ|=9.656 MVar
- transformer active losses: series=0.282 MW, no-load G/shunt=1.54 MW, total=1.822 MW; total active network loss=10.355 MW
- note: ALPHA S1 dV=-0.598 kV; BETA1 S1 dV=-0.462 kV; DELTA1S1 dV=-0.852 kV; DELTA2S1 dV=0.235 kV

### Top voltage deviations

- DELTA1S1: model 389.448 kV vs FOR002 390.3 kV; dV=-0.852 kV; dVa=-0.022 deg
- ALPHA S1: model 389.102 kV vs FOR002 389.7 kV; dV=-0.598 kV; dVa=0.016 deg
- BETA1 S1: model 396.738 kV vs FOR002 397.2 kV; dV=-0.462 kV; dVa=0.049 deg
- DELTA2S1: model 240.635 kV vs FOR002 240.4 kV; dV=0.235 kV; dVa=-0.029 deg
- BETA2 S1: model 231.245 kV vs FOR002 231.4 kV; dV=-0.155 kV; dVa=0.012 deg
- BURG  S1: model 228.511 kV vs FOR002 228.6 kV; dV=-0.089 kV; dVa=0.016 deg
- NORD  S1: model 233.365 kV vs FOR002 233.3 kV; dV=0.065 kV; dVa=0.016 deg
- GAMMA S1: model 226.565 kV vs FOR002 226.6 kV; dV=-0.035 kV; dVa=0.025 deg
- OST   S1: model 229.767 kV vs FOR002 229.8 kV; dV=-0.033 kV; dVa=0.025 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg

### Top branch P deviations

- L1 DELTA2S1 -> WEILERS1: model -73.39 MW vs FOR002 -71.9 MW; dP=-1.49 MW
- L1 BETA2 S1 -> WEILERS1: model -13.762 MW vs FOR002 -15.2 MW; dP=1.438 MW
- T1B DELTA1S1 -> DELTA2S1: model -35.398 MW vs FOR002 -34.4 MW; dP=-0.998 MW
- T1A DELTA1S1 -> DELTA2S1: model -35.398 MW vs FOR002 -34.4 MW; dP=-0.998 MW
- L1A BETA1 S1 -> DELTA1S1: model 26.712 MW vs FOR002 27.4 MW; dP=-0.688 MW
- L1B BETA1 S1 -> DELTA1S1: model 26.224 MW vs FOR002 26.9 MW; dP=-0.676 MW
- T1C BETA1 S1 -> BETA2 S1: model -110.622 MW vs FOR002 -111.0 MW; dP=0.378 MW
- T1A BETA1 S1 -> BETA2 S1: model -110.622 MW vs FOR002 -111.0 MW; dP=0.378 MW
- T1B BETA1 S1 -> BETA2 S1: model -110.622 MW vs FOR002 -111.0 MW; dP=0.378 MW
- L1 ALPHA S1 -> BETA1 S1: model -277.248 MW vs FOR002 -277.6 MW; dP=0.352 MW

### Top branch Q deviations

- L1 DELTA2S1 -> WEILERS1: model 243.456 MVar vs FOR002 233.8 MVar; dQ=9.656 MVar
- L1 BETA2 S1 -> WEILERS1: model -168.082 MVar vs FOR002 -159.7 MVar; dQ=-8.382 MVar
- T1B DELTA1S1 -> DELTA2S1: model 152.738 MVar vs FOR002 147.3 MVar; dQ=5.438 MVar
- T1A DELTA1S1 -> DELTA2S1: model 152.738 MVar vs FOR002 147.3 MVar; dQ=5.438 MVar
- L1A BETA1 S1 -> DELTA1S1: model 55.775 MVar vs FOR002 51.2 MVar; dQ=4.575 MVar
- L1B BETA1 S1 -> DELTA1S1: model 53.541 MVar vs FOR002 49.0 MVar; dQ=4.541 MVar
- L1 BSTADTS1 -> BURG  S1: model 41.407 MVar vs FOR002 37.3 MVar; dQ=4.107 MVar
- T1C BETA1 S1 -> BETA2 S1: model -83.536 MVar vs FOR002 -79.5 MVar; dQ=-4.036 MVar
- T1A BETA1 S1 -> BETA2 S1: model -83.536 MVar vs FOR002 -79.5 MVar; dQ=-4.036 MVar
- T1B BETA1 S1 -> BETA2 S1: model -83.536 MVar vs FOR002 -79.5 MVar; dQ=-4.036 MVar

### Mode winding_over_network

- converged: `true`, iterations: `5`, final mismatch: `1.5422098172789457e-13`
- slack: `BSTADTS1` Sparlectra P/Q=311.876/-218.708; FOR002 P/Q=311.6/-223.2; delta=0.276/4.492
- losses Sparlectra P/Q=10.343/-258.708; FOR002 P/Q=11.6/-263.2; delta=-1.257/4.492
- max deviations: |dV|=2.671 kV (0.007 pu), |dVa|=0.051 deg, branch |dP|=1.741 MW, branch |dQ|=8.868 MVar
- transformer active losses: series=0.282 MW, no-load G/shunt=1.533 MW, total=1.815 MW; total active network loss=10.343 MW
- note: ALPHA S1 dV=-2.441 kV; BETA1 S1 dV=-2.257 kV; DELTA1S1 dV=-2.671 kV; DELTA2S1 dV=0.175 kV

### Top voltage deviations

- DELTA1S1: model 387.629 kV vs FOR002 390.3 kV; dV=-2.671 kV; dVa=-0.021 deg
- ALPHA S1: model 387.259 kV vs FOR002 389.7 kV; dV=-2.441 kV; dVa=0.004 deg
- BETA1 S1: model 394.943 kV vs FOR002 397.2 kV; dV=-2.257 kV; dVa=0.051 deg
- BETA2 S1: model 231.213 kV vs FOR002 231.4 kV; dV=-0.187 kV; dVa=0.014 deg
- DELTA2S1: model 240.575 kV vs FOR002 240.4 kV; dV=0.175 kV; dVa=-0.027 deg
- BURG  S1: model 228.493 kV vs FOR002 228.6 kV; dV=-0.107 kV; dVa=0.017 deg
- GAMMA S1: model 226.547 kV vs FOR002 226.6 kV; dV=-0.053 kV; dVa=0.026 deg
- NORD  S1: model 233.35 kV vs FOR002 233.3 kV; dV=0.05 kV; dVa=0.017 deg
- OST   S1: model 229.76 kV vs FOR002 229.8 kV; dV=-0.04 kV; dVa=0.025 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg

### Top branch P deviations

- L1 DELTA2S1 -> WEILERS1: model -73.513 MW vs FOR002 -71.9 MW; dP=-1.613 MW
- L1 BETA2 S1 -> WEILERS1: model -13.66 MW vs FOR002 -15.2 MW; dP=1.54 MW
- T1B DELTA1S1 -> DELTA2S1: model -35.478 MW vs FOR002 -34.4 MW; dP=-1.078 MW
- T1A DELTA1S1 -> DELTA2S1: model -35.478 MW vs FOR002 -34.4 MW; dP=-1.078 MW
- L1A BETA1 S1 -> DELTA1S1: model 26.649 MW vs FOR002 27.4 MW; dP=-0.751 MW
- L1B BETA1 S1 -> DELTA1S1: model 26.161 MW vs FOR002 26.9 MW; dP=-0.739 MW
- T1C BETA1 S1 -> BETA2 S1: model -110.572 MW vs FOR002 -111.0 MW; dP=0.428 MW
- T1A BETA1 S1 -> BETA2 S1: model -110.572 MW vs FOR002 -111.0 MW; dP=0.428 MW
- T1B BETA1 S1 -> BETA2 S1: model -110.572 MW vs FOR002 -111.0 MW; dP=0.428 MW
- L1 ALPHA S1 -> BETA1 S1: model -277.218 MW vs FOR002 -277.6 MW; dP=0.382 MW

### Top branch Q deviations

- L1 BETA2 S1 -> WEILERS1: model -168.262 MVar vs FOR002 -159.7 MVar; dQ=-8.562 MVar
- L1 DELTA2S1 -> WEILERS1: model 242.218 MVar vs FOR002 233.8 MVar; dQ=8.418 MVar
- L1 BSTADTS1 -> BURG  S1: model 42.306 MVar vs FOR002 37.3 MVar; dQ=5.006 MVar
- L1A BETA1 S1 -> DELTA1S1: model 56.019 MVar vs FOR002 51.2 MVar; dQ=4.819 MVar
- L1B BETA1 S1 -> DELTA1S1: model 53.791 MVar vs FOR002 49.0 MVar; dQ=4.791 MVar
- T1B DELTA1S1 -> DELTA2S1: model 151.93 MVar vs FOR002 147.3 MVar; dQ=4.63 MVar
- T1A DELTA1S1 -> DELTA2S1: model 151.93 MVar vs FOR002 147.3 MVar; dQ=4.63 MVar
- T1C BETA1 S1 -> BETA2 S1: model -83.888 MVar vs FOR002 -79.5 MVar; dQ=-4.388 MVar
- T1A BETA1 S1 -> BETA2 S1: model -83.888 MVar vs FOR002 -79.5 MVar; dQ=-4.388 MVar
- T1B BETA1 S1 -> BETA2 S1: model -83.888 MVar vs FOR002 -79.5 MVar; dQ=-4.388 MVar

## Case D: slack moved to ALPHA

### Mode neutral_one

- converged: `true`, iterations: `5`, final mismatch: `2.0938806244430452e-13`
- slack: `ALPHA S1` Sparlectra P/Q=9.558/-367.918; FOR002 P/Q=9.6/-367.9; delta=-0.042/-0.018
- losses Sparlectra P/Q=7.901/-307.918; FOR002 P/Q=9.6/-307.9; delta=-1.699/-0.018
- max deviations: |dV|=0.034 kV (0.0 pu), |dVa|=0.047 deg, branch |dP|=0.213 MW, branch |dQ|=0.05 MVar
- transformer active losses: series=0.157 MW, no-load G/shunt=1.657 MW, total=1.813 MW; total active network loss=7.901 MW
- note: slack=ALPHA S1; BSTADTS1 P/Q=300.0/100.0.

### Top voltage deviations

- NORD  S1: model 242.534 kV vs FOR002 242.5 kV; dV=0.034 kV; dVa=0.012 deg
- GAMMA S1: model 237.133 kV vs FOR002 237.1 kV; dV=0.033 kV; dVa=-0.004 deg
- DELTA1S1: model 414.474 kV vs FOR002 414.5 kV; dV=-0.026 kV; dVa=-0.016 deg
- ASTADTS1: model 240.924 kV vs FOR002 240.9 kV; dV=0.024 kV; dVa=0.046 deg
- BSTADTS1: model 241.578 kV vs FOR002 241.6 kV; dV=-0.022 kV; dVa=-0.007 deg
- OST   S1: model 240.22 kV vs FOR002 240.2 kV; dV=0.02 kV; dVa=0.038 deg
- BETA2 S1: model 240.684 kV vs FOR002 240.7 kV; dV=-0.016 kV; dVa=0.042 deg
- WEILERS1: model 241.886 kV vs FOR002 241.9 kV; dV=-0.014 kV; dVa=-0.028 deg
- BURG  S1: model 238.989 kV vs FOR002 239.0 kV; dV=-0.011 kV; dVa=0.019 deg
- SUED  S1: model 242.508 kV vs FOR002 242.5 kV; dV=0.008 kV; dVa=0.041 deg

### Top branch P deviations

- T1B DELTA1S1 -> DELTA2S1: model -48.383 MW vs FOR002 -48.2 MW; dP=-0.183 MW
- T1A DELTA1S1 -> DELTA2S1: model -48.383 MW vs FOR002 -48.2 MW; dP=-0.183 MW
- T1C BETA1 S1 -> BETA2 S1: model -99.139 MW vs FOR002 -99.0 MW; dP=-0.139 MW
- T1A BETA1 S1 -> BETA2 S1: model -99.139 MW vs FOR002 -99.0 MW; dP=-0.139 MW
- T1B BETA1 S1 -> BETA2 S1: model -99.139 MW vs FOR002 -99.0 MW; dP=-0.139 MW
- L1 BSTADTS1 -> WEILERS1: model -4.049 MW vs FOR002 -4.0 MW; dP=-0.049 MW
- L1 DELTA2S1 -> ASTADTS1: model -3.951 MW vs FOR002 -4.0 MW; dP=0.049 MW
- L1 ASTADTS1 -> NORD  S1: model -99.746 MW vs FOR002 -99.7 MW; dP=-0.046 MW
- L1 ALPHA S1 -> BETA1 S1: model -263.554 MW vs FOR002 -263.6 MW; dP=0.046 MW
- L1B BETA1 S1 -> DELTA1S1: model 15.545 MW vs FOR002 15.5 MW; dP=0.045 MW

### Top branch Q deviations

- L1 BETA2 S1 -> WEILERS1: model -65.25 MVar vs FOR002 -65.3 MVar; dQ=0.05 MVar
- L1A BETA1 S1 -> DELTA1S1: model -50.046 MVar vs FOR002 -50.0 MVar; dQ=-0.046 MVar
- L1 ASTADTS1 -> OST   S1: model 16.661 MVar vs FOR002 16.7 MVar; dQ=-0.039 MVar
- L1 ALPHA S1 -> DELTA1S1: model -224.13 MVar vs FOR002 -224.1 MVar; dQ=-0.03 MVar
- L1 DELTA2S1 -> ASTADTS1: model -15.329 MVar vs FOR002 -15.3 MVar; dQ=-0.029 MVar
- T1B DELTA1S1 -> DELTA2S1: model -39.576 MVar vs FOR002 -39.6 MVar; dQ=0.024 MVar
- T1A DELTA1S1 -> DELTA2S1: model -39.576 MVar vs FOR002 -39.6 MVar; dQ=0.024 MVar
- L1 BSTADTS1 -> BURG  S1: model 104.378 MVar vs FOR002 104.4 MVar; dQ=-0.022 MVar
- L1 BSTADTS1 -> OST   S1: model 31.022 MVar vs FOR002 31.0 MVar; dQ=0.022 MVar
- L1 BSTADTS1 -> WEILERS1: model -6.978 MVar vs FOR002 -7.0 MVar; dQ=0.022 MVar

### Mode winding_over_network

- converged: `true`, iterations: `5`, final mismatch: `1.3539169785303784e-13`
- slack: `ALPHA S1` Sparlectra P/Q=9.528/-368.706; FOR002 P/Q=9.6/-367.9; delta=-0.072/-0.806
- losses Sparlectra P/Q=7.864/-308.706; FOR002 P/Q=9.6/-307.9; delta=-1.736/-0.806
- max deviations: |dV|=1.095 kV (0.005 pu), |dVa|=0.047 deg, branch |dP|=0.23 MW, branch |dQ|=0.558 MVar
- transformer active losses: series=0.156 MW, no-load G/shunt=1.664 MW, total=1.82 MW; total active network loss=7.864 MW
- note: slack=ALPHA S1; BSTADTS1 P/Q=300.0/100.0.

### Top voltage deviations

- GAMMA S1: model 238.195 kV vs FOR002 237.1 kV; dV=1.095 kV; dVa=0.003 deg
- NORD  S1: model 243.574 kV vs FOR002 242.5 kV; dV=1.074 kV; dVa=0.004 deg
- ASTADTS1: model 241.971 kV vs FOR002 240.9 kV; dV=1.071 kV; dVa=0.043 deg
- OST   S1: model 241.27 kV vs FOR002 240.2 kV; dV=1.07 kV; dVa=0.037 deg
- DELTA2S1: model 240.748 kV vs FOR002 239.7 kV; dV=1.048 kV; dVa=-0.031 deg
- SUED  S1: model 243.548 kV vs FOR002 242.5 kV; dV=1.048 kV; dVa=0.031 deg
- BURG  S1: model 240.043 kV vs FOR002 239.0 kV; dV=1.043 kV; dVa=0.02 deg
- WEILERS1: model 242.926 kV vs FOR002 241.9 kV; dV=1.026 kV; dVa=-0.035 deg
- BETA2 S1: model 241.726 kV vs FOR002 240.7 kV; dV=1.026 kV; dVa=0.035 deg
- BSTADTS1: model 242.623 kV vs FOR002 241.6 kV; dV=1.023 kV; dVa=-0.014 deg

### Top branch P deviations

- T1B DELTA1S1 -> DELTA2S1: model -48.425 MW vs FOR002 -48.2 MW; dP=-0.225 MW
- T1A DELTA1S1 -> DELTA2S1: model -48.425 MW vs FOR002 -48.2 MW; dP=-0.225 MW
- T1C BETA1 S1 -> BETA2 S1: model -99.123 MW vs FOR002 -99.0 MW; dP=-0.123 MW
- T1A BETA1 S1 -> BETA2 S1: model -99.123 MW vs FOR002 -99.0 MW; dP=-0.123 MW
- T1B BETA1 S1 -> BETA2 S1: model -99.123 MW vs FOR002 -99.0 MW; dP=-0.123 MW
- L1 DELTA2S1 -> WEILERS1: model -93.194 MW vs FOR002 -93.1 MW; dP=-0.094 MW
- L1 WEILERS1 -> ASTADTS1: model 158.341 MW vs FOR002 158.4 MW; dP=-0.059 MW
- L1 ASTADTS1 -> NORD  S1: model -99.748 MW vs FOR002 -99.7 MW; dP=-0.048 MW
- L1C BETA2 S1 -> BURG  S1: model 98.655 MW vs FOR002 98.7 MW; dP=-0.045 MW
- L1B BETA2 S1 -> BURG  S1: model 98.655 MW vs FOR002 98.7 MW; dP=-0.045 MW

### Top branch Q deviations

- L1 ALPHA S1 -> BETA1 S1: model -344.306 MVar vs FOR002 -343.8 MVar; dQ=-0.506 MVar
- L1 ALPHA S1 -> DELTA1S1: model -224.4 MVar vs FOR002 -224.1 MVar; dQ=-0.3 MVar
- L1 WEILERS1 -> ASTADTS1: model 46.812 MVar vs FOR002 47.0 MVar; dQ=-0.188 MVar
- T1C BETA1 S1 -> BETA2 S1: model -72.18 MVar vs FOR002 -72.0 MVar; dQ=-0.18 MVar
- T1A BETA1 S1 -> BETA2 S1: model -72.18 MVar vs FOR002 -72.0 MVar; dQ=-0.18 MVar
- T1B BETA1 S1 -> BETA2 S1: model -72.18 MVar vs FOR002 -72.0 MVar; dQ=-0.18 MVar
- L1 DELTA2S1 -> WEILERS1: model -66.256 MVar vs FOR002 -66.1 MVar; dQ=-0.156 MVar
- L1A BETA2 S1 -> BURG  S1: model 13.657 MVar vs FOR002 13.8 MVar; dQ=-0.143 MVar
- L1 BETA2 S1 -> WEILERS1: model -65.436 MVar vs FOR002 -65.3 MVar; dQ=-0.136 MVar
- L1 DELTA2S1 -> ASTADTS1: model -15.434 MVar vs FOR002 -15.3 MVar; dQ=-0.134 MVar

## Case E: DELTA tap step 7 plus 60 degree skew-angle regulator

### Mode neutral_one

- converged: `true`, iterations: `5`, final mismatch: `4.0027129520637083e-13`
- slack: `BSTADTS1` Sparlectra P/Q=311.125/-230.645; FOR002 P/Q=311.0/-231.5; delta=0.125/0.855
- losses Sparlectra P/Q=9.561/-270.645; FOR002 P/Q=11.0/-271.5; delta=-1.439/0.855
- max deviations: |dV|=0.247 kV (0.001 pu), |dVa|=0.036 deg, branch |dP|=2.697 MW, branch |dQ|=2.438 MVar
- transformer active losses: series=0.386 MW, no-load G/shunt=1.564 MW, total=1.95 MW; total active network loss=9.561 MW
- note: BETA transformer P(from) MW=-182.573, -182.573, -182.573; DELTA transformer P(from) MW=72.441, 72.441.

### Top voltage deviations

- DELTA1S1: model 398.353 kV vs FOR002 398.6 kV; dV=-0.247 kV; dVa=0.019 deg
- ALPHA S1: model 395.732 kV vs FOR002 395.9 kV; dV=-0.168 kV; dVa=0.006 deg
- BETA1 S1: model 401.955 kV vs FOR002 402.1 kV; dV=-0.145 kV; dVa=0.021 deg
- GAMMA S1: model 227.326 kV vs FOR002 227.4 kV; dV=-0.074 kV; dVa=0.001 deg
- OST   S1: model 229.655 kV vs FOR002 229.7 kV; dV=-0.045 kV; dVa=0.005 deg
- NORD  S1: model 233.141 kV vs FOR002 233.1 kV; dV=0.041 kV; dVa=-0.023 deg
- BURG  S1: model 229.265 kV vs FOR002 229.3 kV; dV=-0.035 kV; dVa=-0.013 deg
- DELTA2S1: model 238.232 kV vs FOR002 238.2 kV; dV=0.032 kV; dVa=0.019 deg
- ASTADTS1: model 231.472 kV vs FOR002 231.5 kV; dV=-0.028 kV; dVa=-0.036 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg

### Top branch P deviations

- L1 DELTA2S1 -> WEILERS1: model 108.797 MW vs FOR002 106.1 MW; dP=2.697 MW
- L1 BETA2 S1 -> WEILERS1: model -170.247 MW vs FOR002 -168.0 MW; dP=-2.247 MW
- T1B DELTA1S1 -> DELTA2S1: model 72.441 MW vs FOR002 70.9 MW; dP=1.541 MW
- T1A DELTA1S1 -> DELTA2S1: model 72.441 MW vs FOR002 70.9 MW; dP=1.541 MW
- T1C BETA1 S1 -> BETA2 S1: model -182.573 MW vs FOR002 -181.3 MW; dP=-1.273 MW
- T1A BETA1 S1 -> BETA2 S1: model -182.573 MW vs FOR002 -181.3 MW; dP=-1.273 MW
- T1B BETA1 S1 -> BETA2 S1: model -182.573 MW vs FOR002 -181.3 MW; dP=-1.273 MW
- L1B BETA1 S1 -> DELTA1S1: model 107.033 MW vs FOR002 105.8 MW; dP=1.233 MW
- L1A BETA1 S1 -> DELTA1S1: model 108.997 MW vs FOR002 107.8 MW; dP=1.197 MW
- L1 BSTADTS1 -> BURG  S1: model 281.701 MW vs FOR002 280.7 MW; dP=1.001 MW

### Top branch Q deviations

- L1 DELTA2S1 -> WEILERS1: model 124.238 MVar vs FOR002 121.8 MVar; dQ=2.438 MVar
- L1 BETA2 S1 -> WEILERS1: model -66.46 MVar vs FOR002 -64.2 MVar; dQ=-2.26 MVar
- T1B DELTA1S1 -> DELTA2S1: model 79.84 MVar vs FOR002 78.5 MVar; dQ=1.34 MVar
- T1A DELTA1S1 -> DELTA2S1: model 79.84 MVar vs FOR002 78.5 MVar; dQ=1.34 MVar
- L1A BETA1 S1 -> DELTA1S1: model -2.893 MVar vs FOR002 -4.1 MVar; dQ=1.207 MVar
- L1B BETA1 S1 -> DELTA1S1: model -4.098 MVar vs FOR002 -5.3 MVar; dQ=1.202 MVar
- L1 BSTADTS1 -> BURG  S1: model -3.739 MVar vs FOR002 -4.9 MVar; dQ=1.161 MVar
- T1C BETA1 S1 -> BETA2 S1: model -32.02 MVar vs FOR002 -30.9 MVar; dQ=-1.12 MVar
- T1A BETA1 S1 -> BETA2 S1: model -32.02 MVar vs FOR002 -30.9 MVar; dQ=-1.12 MVar
- T1B BETA1 S1 -> BETA2 S1: model -32.02 MVar vs FOR002 -30.9 MVar; dQ=-1.12 MVar

### Mode winding_over_network

- converged: `true`, iterations: `5`, final mismatch: `2.868554827115666e-13`
- slack: `BSTADTS1` Sparlectra P/Q=311.107/-228.056; FOR002 P/Q=311.0/-231.5; delta=0.107/3.444
- losses Sparlectra P/Q=9.55/-268.056; FOR002 P/Q=11.0/-271.5; delta=-1.45/3.444
- max deviations: |dV|=2.095 kV (0.005 pu), |dVa|=0.036 deg, branch |dP|=2.176 MW, branch |dQ|=2.805 MVar
- transformer active losses: series=0.385 MW, no-load G/shunt=1.557 MW, total=1.942 MW; total active network loss=9.55 MW
- note: BETA transformer P(from) MW=-182.368, -182.368, -182.368; DELTA transformer P(from) MW=72.128, 72.128.

### Top voltage deviations

- DELTA1S1: model 396.505 kV vs FOR002 398.6 kV; dV=-2.095 kV; dVa=0.016 deg
- ALPHA S1: model 393.855 kV vs FOR002 395.9 kV; dV=-2.045 kV; dVa=-0.005 deg
- BETA1 S1: model 400.123 kV vs FOR002 402.1 kV; dV=-1.977 kV; dVa=0.026 deg
- GAMMA S1: model 227.305 kV vs FOR002 227.4 kV; dV=-0.095 kV; dVa=0.002 deg
- BURG  S1: model 229.245 kV vs FOR002 229.3 kV; dV=-0.055 kV; dVa=-0.011 deg
- OST   S1: model 229.648 kV vs FOR002 229.7 kV; dV=-0.052 kV; dVa=0.005 deg
- ASTADTS1: model 231.458 kV vs FOR002 231.5 kV; dV=-0.042 kV; dVa=-0.036 deg
- BETA2 S1: model 232.465 kV vs FOR002 232.5 kV; dV=-0.035 kV; dVa=0.005 deg
- SUED  S1: model 230.973 kV vs FOR002 231.0 kV; dV=-0.027 kV; dVa=-0.016 deg
- NORD  S1: model 233.127 kV vs FOR002 233.1 kV; dV=0.027 kV; dVa=-0.023 deg

### Top branch P deviations

- L1 DELTA2S1 -> WEILERS1: model 108.276 MW vs FOR002 106.1 MW; dP=2.176 MW
- L1 BETA2 S1 -> WEILERS1: model -169.799 MW vs FOR002 -168.0 MW; dP=-1.799 MW
- T1B DELTA1S1 -> DELTA2S1: model 72.128 MW vs FOR002 70.9 MW; dP=1.228 MW
- T1A DELTA1S1 -> DELTA2S1: model 72.128 MW vs FOR002 70.9 MW; dP=1.228 MW
- T1C BETA1 S1 -> BETA2 S1: model -182.368 MW vs FOR002 -181.3 MW; dP=-1.068 MW
- T1A BETA1 S1 -> BETA2 S1: model -182.368 MW vs FOR002 -181.3 MW; dP=-1.068 MW
- T1B BETA1 S1 -> BETA2 S1: model -182.368 MW vs FOR002 -181.3 MW; dP=-1.068 MW
- L1B BETA1 S1 -> DELTA1S1: model 106.801 MW vs FOR002 105.8 MW; dP=1.001 MW
- L1A BETA1 S1 -> DELTA1S1: model 108.761 MW vs FOR002 107.8 MW; dP=0.961 MW
- L1 BSTADTS1 -> BURG  S1: model 281.528 MW vs FOR002 280.7 MW; dP=0.828 MW

### Top branch Q deviations

- L1 BETA2 S1 -> WEILERS1: model -66.903 MVar vs FOR002 -64.2 MVar; dQ=-2.703 MVar
- L1 BSTADTS1 -> BURG  S1: model -2.716 MVar vs FOR002 -4.9 MVar; dQ=2.184 MVar
- T1C BETA1 S1 -> BETA2 S1: model -32.52 MVar vs FOR002 -30.9 MVar; dQ=-1.62 MVar
- T1A BETA1 S1 -> BETA2 S1: model -32.52 MVar vs FOR002 -30.9 MVar; dQ=-1.62 MVar
- T1B BETA1 S1 -> BETA2 S1: model -32.52 MVar vs FOR002 -30.9 MVar; dQ=-1.62 MVar
- L1A BETA1 S1 -> DELTA1S1: model -2.487 MVar vs FOR002 -4.1 MVar; dQ=1.613 MVar
- L1B BETA1 S1 -> DELTA1S1: model -3.688 MVar vs FOR002 -5.3 MVar; dQ=1.612 MVar
- L1 DELTA2S1 -> WEILERS1: model 123.277 MVar vs FOR002 121.8 MVar; dQ=1.477 MVar
- L1 WEILERS1 -> ASTADTS1: model 191.557 MVar vs FOR002 192.7 MVar; dQ=-1.143 MVar
- L1 ALPHA S1 -> DELTA1S1: model -66.387 MVar vs FOR002 -67.3 MVar; dQ=0.913 MVar

