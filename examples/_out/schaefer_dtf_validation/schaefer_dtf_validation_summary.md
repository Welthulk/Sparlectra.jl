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

- converged: `false`, iterations: `1`, final mismatch: `NaN`
- slack: `BSTADTS1` Sparlectra P/Q=0.0/-3.11098874593e8; FOR002 P/Q=309.2/-201.5; delta=-309.2/-3.11098673093e8
- losses Sparlectra P/Q=n/a/n/a; FOR002 P/Q=9.2/-308.0; delta=n/a/n/a
- max deviations: |dV|=4.3 kV (107.5 pu), |dVa|=1.9 deg, branch |dP|=1.16918355759e8 MW, branch |dQ|=6.46203827762e8 MVar
- transformer active losses: series=n/a MW, no-load G/shunt=n/a MW, total=n/a MW; total active network loss=n/a MW
- note: NORD effective type=PV, V=230.5 kV, Q=5.14123624341e8 MVar.

### Top voltage deviations

- BETA2 S1: model 230.0 kV vs FOR002 234.3 kV; dV=-4.3 kV; dVa=0.1 deg
- DELTA2S1: model 230.0 kV vs FOR002 234.0 kV; dV=-4.0 kV; dVa=0.8 deg
- WEILERS1: model 230.0 kV vs FOR002 233.6 kV; dV=-3.6 kV; dVa=0.1 deg
- GAMMA S1: model 230.0 kV vs FOR002 228.3 kV; dV=1.7 kV; dVa=1.9 deg
- SUED  S1: model 230.0 kV vs FOR002 231.0 kV; dV=-1.0 kV; dVa=-0.4 deg
- OST   S1: model 230.0 kV vs FOR002 229.2 kV; dV=0.8 kV; dVa=0.8 deg
- ASTADTS1: model 230.0 kV vs FOR002 230.6 kV; dV=-0.6 kV; dVa=0.5 deg
- BURG  S1: model 230.0 kV vs FOR002 230.2 kV; dV=-0.2 kV; dVa=1.1 deg
- NORD  S1: model 230.5 kV vs FOR002 230.5 kV; dV=0.0 kV; dVa=-0.2 deg
- BSTADTS1: model 230.0 kV vs FOR002 230.0 kV; dV=0.0 kV; dVa=0.0 deg

### Top branch P deviations

- L1 ASTADTS1 -> NORD  S1: model -1.16664836549e8 MW vs FOR002 -99.7 MW; dP=-1.16664736849e8 MW
- L1 BSTADTS1 -> BURG  S1: model 0.0 MW vs FOR002 210.8 MW; dP=-210.8 MW
- L1 WEILERS1 -> ASTADTS1: model 0.0 MW vs FOR002 153.4 MW; dP=-153.4 MW
- L1B BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.6 MW; dP=-150.6 MW
- L1A BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.6 MW; dP=-150.6 MW
- L1 BSTADTS1 -> OST   S1: model 0.0 MW vs FOR002 136.2 MW; dP=-136.2 MW
- L1C BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 97.8 MW; dP=-97.8 MW
- L1B BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 97.8 MW; dP=-97.8 MW
- L1A BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 97.4 MW; dP=-97.4 MW
- L1 DELTA2S1 -> WEILERS1: model 0.0 MW vs FOR002 -93.6 MW; dP=93.6 MW

### Top branch Q deviations

- L1 ASTADTS1 -> NORD  S1: model -6.46203804062e8 MVar vs FOR002 23.7 MVar; dQ=-6.46203827762e8 MVar
- L1 DELTA2S1 -> ASTADTS1: model -2.35574020182e8 MVar vs FOR002 16.4 MVar; dQ=-2.35574036582e8 MVar
- L1 BSTADTS1 -> WEILERS1: model -1.41167952169e8 MVar vs FOR002 -39.8 MVar; dQ=-1.41167912369e8 MVar
- L1 DELTA2S1 -> WEILERS1: model -6.652539746e7 MVar vs FOR002 31.3 MVar; dQ=-6.652542876e7 MVar
- L1A BETA2 S1 -> BURG  S1: model -6.5731327729e7 MVar vs FOR002 64.5 MVar; dQ=-6.5731392229e7 MVar
- L1C BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 63.5 MVar; dQ=-6.5378471348e7 MVar
- L1B BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 63.5 MVar; dQ=-6.5378471348e7 MVar
- L1 BETA2 S1 -> WEILERS1: model -4.8526483558e7 MVar vs FOR002 32.9 MVar; dQ=-4.8526516458e7 MVar
- L1 BSTADTS1 -> BURG  S1: model -4.8526483558e7 MVar vs FOR002 -43.5 MVar; dQ=-4.8526440058e7 MVar
- L1B BURG  S1 -> GAMMA S1: model -3.352738864e7 MVar vs FOR002 50.5 MVar; dQ=-3.352743914e7 MVar

## Case C: DELTA fixed tap step 7

- converged: `false`, iterations: `1`, final mismatch: `NaN`
- slack: `BSTADTS1` Sparlectra P/Q=0.0/-3.11098874593e8; FOR002 P/Q=311.6/-223.2; delta=-311.6/-3.11098651393e8
- losses Sparlectra P/Q=n/a/n/a; FOR002 P/Q=11.6/-263.2; delta=n/a/n/a
- max deviations: |dV|=10.4 kV (260.0 pu), |dVa|=1.9 deg, branch |dP|=222.4 MW, branch |dQ|=2.35574072682e8 MVar
- transformer active losses: series=n/a MW, no-load G/shunt=n/a MW, total=n/a MW; total active network loss=n/a MW
- note: ALPHA S1 dV=n/a kV; BETA1 S1 dV=n/a kV; DELTA1S1 dV=n/a kV; DELTA2S1 dV=-10.4 kV

### Top voltage deviations

- DELTA2S1: model 230.0 kV vs FOR002 240.4 kV; dV=-10.4 kV; dVa=0.9 deg
- WEILERS1: model 230.0 kV vs FOR002 234.6 kV; dV=-4.6 kV; dVa=0.1 deg
- GAMMA S1: model 230.0 kV vs FOR002 226.6 kV; dV=3.4 kV; dVa=1.9 deg
- NORD  S1: model 230.0 kV vs FOR002 233.3 kV; dV=-3.3 kV; dVa=-0.1 deg
- ASTADTS1: model 230.0 kV vs FOR002 231.7 kV; dV=-1.7 kV; dVa=0.5 deg
- BETA2 S1: model 230.0 kV vs FOR002 231.4 kV; dV=-1.4 kV; dVa=0.1 deg
- BURG  S1: model 230.0 kV vs FOR002 228.6 kV; dV=1.4 kV; dVa=1.1 deg
- SUED  S1: model 230.0 kV vs FOR002 231.0 kV; dV=-1.0 kV; dVa=-0.4 deg
- OST   S1: model 230.0 kV vs FOR002 229.8 kV; dV=0.2 kV; dVa=0.8 deg
- BSTADTS1: model 230.0 kV vs FOR002 230.0 kV; dV=0.0 kV; dVa=0.0 deg

### Top branch P deviations

- L1 BSTADTS1 -> BURG  S1: model 0.0 MW vs FOR002 222.4 MW; dP=-222.4 MW
- L1 WEILERS1 -> ASTADTS1: model 0.0 MW vs FOR002 156.2 MW; dP=-156.2 MW
- L1B BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.6 MW; dP=-150.6 MW
- L1A BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.6 MW; dP=-150.6 MW
- L1 BSTADTS1 -> OST   S1: model 0.0 MW vs FOR002 133.5 MW; dP=-133.5 MW
- L1 ASTADTS1 -> NORD  S1: model 0.0 MW vs FOR002 -99.7 MW; dP=99.7 MW
- L1C BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 93.7 MW; dP=-93.7 MW
- L1B BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 93.7 MW; dP=-93.7 MW
- L1A BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 93.5 MW; dP=-93.5 MW
- L1B BSTADTS1 -> SUED  S1: model 0.0 MW vs FOR002 -74.9 MW; dP=74.9 MW

### Top branch Q deviations

- L1 DELTA2S1 -> ASTADTS1: model -2.35574020182e8 MVar vs FOR002 52.5 MVar; dQ=-2.35574072682e8 MVar
- L1 BSTADTS1 -> WEILERS1: model -1.41167952169e8 MVar vs FOR002 -49.0 MVar; dQ=-1.41167903169e8 MVar
- L1 DELTA2S1 -> WEILERS1: model -6.652539746e7 MVar vs FOR002 233.8 MVar; dQ=-6.652563126e7 MVar
- L1 ASTADTS1 -> NORD  S1: model -6.652539746e7 MVar vs FOR002 -42.7 MVar; dQ=-6.652535476e7 MVar
- L1A BETA2 S1 -> BURG  S1: model -6.5731327729e7 MVar vs FOR002 37.1 MVar; dQ=-6.5731364829e7 MVar
- L1C BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 36.2 MVar; dQ=-6.5378444048e7 MVar
- L1B BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 36.2 MVar; dQ=-6.5378444048e7 MVar
- L1 BSTADTS1 -> BURG  S1: model -4.8526483558e7 MVar vs FOR002 37.3 MVar; dQ=-4.8526520858e7 MVar
- L1 BETA2 S1 -> WEILERS1: model -4.8526483558e7 MVar vs FOR002 -159.7 MVar; dQ=-4.8526323858e7 MVar
- L1 ASTADTS1 -> OST   S1: model -3.352738864e7 MVar vs FOR002 69.4 MVar; dQ=-3.352745804e7 MVar

## Case D: slack moved to ALPHA

- converged: `false`, iterations: `1`, final mismatch: `NaN`
- slack: `ALPHA S1` Sparlectra P/Q=n/a/n/a; FOR002 P/Q=9.6/-367.9; delta=n/a/n/a
- losses Sparlectra P/Q=n/a/n/a; FOR002 P/Q=9.6/-307.9; delta=n/a/n/a
- max deviations: |dV|=12.5 kV (312.5 pu), |dVa|=2.2 deg, branch |dP|=207.0 MW, branch |dQ|=2.35574019982e8 MVar
- transformer active losses: series=n/a MW, no-load G/shunt=n/a MW, total=n/a MW; total active network loss=n/a MW
- note: slack=ALPHA S1; BSTADTS1 P/Q=300.0/100.0.

### Top voltage deviations

- NORD  S1: model 230.0 kV vs FOR002 242.5 kV; dV=-12.5 kV; dVa=-2.1 deg
- SUED  S1: model 230.0 kV vs FOR002 242.5 kV; dV=-12.5 kV; dVa=-2.2 deg
- WEILERS1: model 230.0 kV vs FOR002 241.9 kV; dV=-11.9 kV; dVa=-2.0 deg
- BSTADTS1: model 230.0 kV vs FOR002 241.6 kV; dV=-11.6 kV; dVa=-1.9 deg
- ASTADTS1: model 230.0 kV vs FOR002 240.9 kV; dV=-10.9 kV; dVa=-1.5 deg
- BETA2 S1: model 230.0 kV vs FOR002 240.7 kV; dV=-10.7 kV; dVa=-2.0 deg
- OST   S1: model 230.0 kV vs FOR002 240.2 kV; dV=-10.2 kV; dVa=-1.2 deg
- DELTA2S1: model 230.0 kV vs FOR002 239.7 kV; dV=-9.7 kV; dVa=-1.5 deg
- BURG  S1: model 230.0 kV vs FOR002 239.0 kV; dV=-9.0 kV; dVa=-1.0 deg
- GAMMA S1: model 230.0 kV vs FOR002 237.1 kV; dV=-7.1 kV; dVa=-0.3 deg

### Top branch P deviations

- L1 BSTADTS1 -> BURG  S1: model 0.0 MW vs FOR002 207.0 MW; dP=-207.0 MW
- L1 WEILERS1 -> ASTADTS1: model 0.0 MW vs FOR002 158.4 MW; dP=-158.4 MW
- L1B BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.5 MW; dP=-150.5 MW
- L1A BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.5 MW; dP=-150.5 MW
- L1 BSTADTS1 -> OST   S1: model 0.0 MW vs FOR002 136.0 MW; dP=-136.0 MW
- L1 ASTADTS1 -> NORD  S1: model 0.0 MW vs FOR002 -99.7 MW; dP=99.7 MW
- L1C BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 98.7 MW; dP=-98.7 MW
- L1A BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 98.7 MW; dP=-98.7 MW
- L1B BETA2 S1 -> BURG  S1: model 0.0 MW vs FOR002 98.7 MW; dP=-98.7 MW
- L1 DELTA2S1 -> WEILERS1: model 0.0 MW vs FOR002 -93.1 MW; dP=93.1 MW

### Top branch Q deviations

- L1 DELTA2S1 -> ASTADTS1: model -2.35574020182e8 MVar vs FOR002 -15.3 MVar; dQ=-2.35574004882e8 MVar
- L1 BSTADTS1 -> WEILERS1: model -1.41167952169e8 MVar vs FOR002 -7.0 MVar; dQ=-1.41167945169e8 MVar
- L1 ASTADTS1 -> NORD  S1: model -6.652539746e7 MVar vs FOR002 -43.2 MVar; dQ=-6.652535426e7 MVar
- L1 DELTA2S1 -> WEILERS1: model -6.652539746e7 MVar vs FOR002 -66.1 MVar; dQ=-6.652533136e7 MVar
- L1A BETA2 S1 -> BURG  S1: model -6.5731327729e7 MVar vs FOR002 13.8 MVar; dQ=-6.5731341529e7 MVar
- L1C BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 12.9 MVar; dQ=-6.5378420748e7 MVar
- L1B BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 12.9 MVar; dQ=-6.5378420748e7 MVar
- L1 BSTADTS1 -> BURG  S1: model -4.8526483558e7 MVar vs FOR002 104.4 MVar; dQ=-4.8526587958e7 MVar
- L1 BETA2 S1 -> WEILERS1: model -4.8526483558e7 MVar vs FOR002 -65.3 MVar; dQ=-4.8526418258e7 MVar
- L1B BURG  S1 -> GAMMA S1: model -3.352738864e7 MVar vs FOR002 50.1 MVar; dQ=-3.352743874e7 MVar

## Case E: DELTA tap step 7 plus 60 degree phase shifter

- converged: `false`, iterations: `1`, final mismatch: `NaN`
- slack: `BSTADTS1` Sparlectra P/Q=0.0/-3.11098874593e8; FOR002 P/Q=311.0/-231.5; delta=-311.0/-3.11098643093e8
- losses Sparlectra P/Q=n/a/n/a; FOR002 P/Q=11.0/-271.5; delta=n/a/n/a
- max deviations: |dV|=8.2 kV (205.0 pu), |dVa|=2.2 deg, branch |dP|=280.7 MW, branch |dQ|=2.35574050782e8 MVar
- transformer active losses: series=n/a MW, no-load G/shunt=n/a MW, total=n/a MW; total active network loss=n/a MW
- note: BETA transformer P(from) MW=n/a, n/a, n/a; DELTA transformer P(from) MW=n/a, n/a.

### Top voltage deviations

- DELTA2S1: model 230.0 kV vs FOR002 238.2 kV; dV=-8.2 kV; dVa=-0.6 deg
- WEILERS1: model 230.0 kV vs FOR002 234.3 kV; dV=-4.3 kV; dVa=-0.1 deg
- NORD  S1: model 230.0 kV vs FOR002 233.1 kV; dV=-3.1 kV; dVa=-0.3 deg
- GAMMA S1: model 230.0 kV vs FOR002 227.4 kV; dV=2.6 kV; dVa=2.2 deg
- BETA2 S1: model 230.0 kV vs FOR002 232.5 kV; dV=-2.5 kV; dVa=0.7 deg
- ASTADTS1: model 230.0 kV vs FOR002 231.5 kV; dV=-1.5 kV; dVa=0.3 deg
- SUED  S1: model 230.0 kV vs FOR002 231.0 kV; dV=-1.0 kV; dVa=-0.4 deg
- BURG  S1: model 230.0 kV vs FOR002 229.3 kV; dV=0.7 kV; dVa=1.4 deg
- OST   S1: model 230.0 kV vs FOR002 229.7 kV; dV=0.3 kV; dVa=0.7 deg
- BSTADTS1: model 230.0 kV vs FOR002 230.0 kV; dV=0.0 kV; dVa=0.0 deg

### Top branch P deviations

- L1 BSTADTS1 -> BURG  S1: model 0.0 MW vs FOR002 280.7 MW; dP=-280.7 MW
- L1 WEILERS1 -> ASTADTS1: model 0.0 MW vs FOR002 174.3 MW; dP=-174.3 MW
- L1 BETA2 S1 -> WEILERS1: model 0.0 MW vs FOR002 -168.0 MW; dP=168.0 MW
- L1B BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.6 MW; dP=-150.6 MW
- L1A BURG  S1 -> GAMMA S1: model 0.0 MW vs FOR002 150.6 MW; dP=-150.6 MW
- L1 BSTADTS1 -> OST   S1: model 0.0 MW vs FOR002 121.4 MW; dP=-121.4 MW
- L1 DELTA2S1 -> WEILERS1: model 0.0 MW vs FOR002 106.1 MW; dP=-106.1 MW
- L1 ASTADTS1 -> NORD  S1: model 0.0 MW vs FOR002 -99.7 MW; dP=99.7 MW
- L1 ASTADTS1 -> OST   S1: model 0.0 MW vs FOR002 79.2 MW; dP=-79.2 MW
- L1B BSTADTS1 -> SUED  S1: model 0.0 MW vs FOR002 -74.9 MW; dP=74.9 MW

### Top branch Q deviations

- L1 DELTA2S1 -> ASTADTS1: model -2.35574020182e8 MVar vs FOR002 30.6 MVar; dQ=-2.35574050782e8 MVar
- L1 BSTADTS1 -> WEILERS1: model -1.41167952169e8 MVar vs FOR002 -44.4 MVar; dQ=-1.41167907769e8 MVar
- L1 DELTA2S1 -> WEILERS1: model -6.652539746e7 MVar vs FOR002 121.8 MVar; dQ=-6.652551926e7 MVar
- L1 ASTADTS1 -> NORD  S1: model -6.652539746e7 MVar vs FOR002 -42.7 MVar; dQ=-6.652535476e7 MVar
- L1A BETA2 S1 -> BURG  S1: model -6.5731327729e7 MVar vs FOR002 51.4 MVar; dQ=-6.5731379129e7 MVar
- L1C BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 50.7 MVar; dQ=-6.5378458548e7 MVar
- L1B BETA2 S1 -> BURG  S1: model -6.5378407848e7 MVar vs FOR002 50.7 MVar; dQ=-6.5378458548e7 MVar
- L1 BSTADTS1 -> BURG  S1: model -4.8526483558e7 MVar vs FOR002 -4.9 MVar; dQ=-4.8526478658e7 MVar
- L1 BETA2 S1 -> WEILERS1: model -4.8526483558e7 MVar vs FOR002 -64.2 MVar; dQ=-4.8526419358e7 MVar
- L1 ASTADTS1 -> OST   S1: model -3.352738864e7 MVar vs FOR002 61.8 MVar; dQ=-3.352745044e7 MVar

