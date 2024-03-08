% This file was generated by Sparlectra exportMatPower.jl
function mpc = bsp3
%% MATPOWER Case Format : Version 2
mpc.version = '2';
%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
mpc.bus = [
%bus	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
1	1	100.0	30.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: ASTADT (ID: 9e0ec64e-6804-42c1-ad1d-dad47dffc3bf kIdx: 1)
2	2	0.0	0.0	0.0	0.0	1	1.03	0.0	110.0	1	1.1	0.9; %Bus: STADION1 (ID: 16a051b1-f70b-4a7c-a185-f2ae05de4695 kIdx: 2)
3	3	0.0	0.0	0.0	0.0	1	1.02	0.0	110.0	1	1.1	0.9; %Bus: VERBUND (slack) (ID: dde3c21f-49ba-4e10-8336-2afd102ee63b kIdx: 3)
];
%% generator data
mpc.gen = [
%bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
2	70.0	0	50.0	-50.0	1.03	100.0	1	0	0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0; %Gen: SynchronousMachine_GEN_STADION1 (ID: 16a051b1-f70b-4a7c-a185-f2ae05de4695)
3	0	0	0	0	1.02	100.0	1	0	0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0; %Gen: ExternalNetworkInjection_EXT_VERBUND (ID: dde3c21f-49ba-4e10-8336-2afd102ee63b)
];
%% branch data
mpc.branch = [
%fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
1	2	0.0	0.08264463	0.00907567	0	0	0	0.0	0.0	1	-360	360; %Branch: Name: L1_2, Typ: LineC (ID: a80b25fc-3186-410e-b3c2-2b7d88367870)
2	3	0.0	0.08264463	0.00907567	0	0	0	0.0	0.0	1	-360	360; %Branch: Name: L2_3, Typ: LineC (ID: 24c88079-1804-4414-8a83-fe57ce3e63a7)
3	1	0.0	0.08264463	0.00907567	0	0	0	0.0	0.0	1	-360	360; %Branch: Name: L3_1, Typ: LineC (ID: 06186ebe-7a9b-4ea2-a741-39c3e2516a12)
];
