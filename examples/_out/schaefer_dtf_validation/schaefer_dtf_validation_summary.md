# Schäfer DTF/FOR002 validation summary

Generated with native `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!`; this is a diagnostic report, not an electrical-model fix.

## Case A: base case

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

## Case E: DELTA tap step 7 plus 60 degree phase shifter

- converged: `false`, iterations: `50`, final mismatch: `23.669309729144167`
- slack: `BSTADTS1` Sparlectra P/Q=1270.097/6962.898; FOR002 P/Q=311.0/-231.5; delta=959.097/7194.398
- losses Sparlectra P/Q=1229.953/9289.658; FOR002 P/Q=11.0/-271.5; delta=1218.953/9561.158
- max deviations: |dV|=398.6 kV (0.996 pu), |dVa|=172.762 deg, branch |dP|=344.499 MW, branch |dQ|=2617.367 MVar
- transformer active losses: series=28.91 MW, no-load G/shunt=0.293 MW, total=29.203 MW; total active network loss=1229.953 MW
- note: BETA transformer P(from) MW=-174.7, -174.7, -174.7; DELTA transformer P(from) MW=0.0, 0.0.

### Top voltage deviations

- DELTA1S1: model 0.0 kV vs FOR002 398.6 kV; dV=-398.6 kV; dVa=-172.762 deg
- BETA1 S1: model 169.822 kV vs FOR002 402.1 kV; dV=-232.278 kV; dVa=-0.786 deg
- ALPHA S1: model 213.129 kV vs FOR002 395.9 kV; dV=-182.771 kV; dVa=-7.956 deg
- DELTA2S1: model 68.072 kV vs FOR002 238.2 kV; dV=-170.128 kV; dVa=7.886 deg
- BETA2 S1: model 138.617 kV vs FOR002 232.5 kV; dV=-93.883 kV; dVa=3.265 deg
- WEILERS1: model 153.424 kV vs FOR002 234.3 kV; dV=-80.876 kV; dVa=2.456 deg
- GAMMA S1: model 173.343 kV vs FOR002 227.4 kV; dV=-54.057 kV; dVa=-0.352 deg
- BURG  S1: model 175.9 kV vs FOR002 229.3 kV; dV=-53.4 kV; dVa=0.184 deg
- ASTADTS1: model 188.716 kV vs FOR002 231.5 kV; dV=-42.784 kV; dVa=0.751 deg
- NORD  S1: model 190.725 kV vs FOR002 233.1 kV; dV=-42.375 kV; dVa=1.071 deg

### Top branch P deviations

- L1 BSTADTS1 -> BURG  S1: model 579.774 MW vs FOR002 280.7 MW; dP=299.074 MW
- L1 ALPHA S1 -> DELTA1S1: model 160.546 MW vs FOR002 -71.0 MW; dP=231.546 MW
- L1 WEILERS1 -> ASTADTS1: model -12.304 MW vs FOR002 174.3 MW; dP=-186.604 MW
- L1 BSTADTS1 -> OST   S1: model 284.519 MW vs FOR002 121.4 MW; dP=163.119 MW
- L1B ASTADTS1 -> BSTADTS1: model -174.904 MW vs FOR002 -35.5 MW; dP=-139.404 MW
- L1A ASTADTS1 -> BSTADTS1: model -174.904 MW vs FOR002 -35.5 MW; dP=-139.404 MW
- L1 ASTADTS1 -> OST   S1: model -46.989 MW vs FOR002 79.2 MW; dP=-126.189 MW
- L1 BSTADTS1 -> WEILERS1: model 107.429 MW vs FOR002 -12.7 MW; dP=120.129 MW
- L1 DELTA2S1 -> WEILERS1: model -12.083 MW vs FOR002 106.1 MW; dP=-118.183 MW
- L1 BETA2 S1 -> WEILERS1: model -64.833 MW vs FOR002 -168.0 MW; dP=103.167 MW

### Top branch Q deviations

- L1 BSTADTS1 -> BURG  S1: model 2612.467 MVar vs FOR002 -4.9 MVar; dQ=2617.367 MVar
- L1 WEILERS1 -> ASTADTS1: model -1848.967 MVar vs FOR002 192.7 MVar; dQ=-2041.667 MVar
- L1 ALPHA S1 -> DELTA1S1: model 1532.255 MVar vs FOR002 -67.3 MVar; dQ=1599.555 MVar
- L1B ASTADTS1 -> BSTADTS1: model -1134.594 MVar vs FOR002 57.9 MVar; dQ=-1192.494 MVar
- L1A ASTADTS1 -> BSTADTS1: model -1134.594 MVar vs FOR002 57.9 MVar; dQ=-1192.494 MVar
- L1 DELTA2S1 -> WEILERS1: model -910.383 MVar vs FOR002 121.8 MVar; dQ=-1032.183 MVar
- L1A BETA1 S1 -> DELTA1S1: model 959.205 MVar vs FOR002 -4.1 MVar; dQ=963.305 MVar
- L1B BETA1 S1 -> DELTA1S1: model 941.733 MVar vs FOR002 -5.3 MVar; dQ=947.033 MVar
- L1 BSTADTS1 -> OST   S1: model 904.717 MVar vs FOR002 -13.4 MVar; dQ=918.117 MVar
- L1 BSTADTS1 -> WEILERS1: model 737.56 MVar vs FOR002 -44.4 MVar; dQ=781.96 MVar

