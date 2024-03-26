% This file was generated by sparlectra exportMatPower.jl
function mpc = bsp2
%% MATPOWER Case Format : Version 2
mpc.version = '2';
%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
mpc.bus = [
%bus	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
1	1	5.0	4.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: K1 (ID: 8fde1ebb-1d8d-4414-b17f-8e6615a747d5 kIdx: 1)
2	1	1.0	1.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: K2 (ID: 005224ad-2fa7-4851-a59a-1a856c716321 kIdx: 2)
3	1	0.0	0.3	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: K3 (ID: 37eaca9e-8cfe-42b0-8357-b4db7b417e3e kIdx: 3)
4	3	0.0	0.0	0.0	0.0	1	1.1	0.0	110.0	1	1.1	0.9; %Bus: K4 (slack) (ID: 6b4c022d-d8db-4bbe-a370-4b6ea05ad210 kIdx: 4)
];
%% generator data
mpc.gen = [
%bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
4	0.0	0	30.0	-30.0	1.1	100.0	1	0	0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0; %Gen: SynchronousMachine_VGEN_#1 (ID: 6b4c022d-d8db-4bbe-a370-4b6ea05ad210)
];
%% branch data
mpc.branch = [
%fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
1	2	0.07136364	0.18636364	0.0	0	0	0	0.0	0.0	1	-360.0	360.0; %Branch: Name: L1_2, Typ: Line (ID: 7cba9ad1-8382-4692-9340-8265dfc98d0c)
1	3	0.03568182	0.09318182	0.0	0	0	0	0.0	0.0	1	-360.0	360.0; %Branch: Name: L1_3, Typ: Line (ID: 9488fa02-a9c4-4f4f-a584-83b23b0adb12)
1	4	0.03568182	0.09318182	0.0	0	0	0	0.0	0.0	1	-360.0	360.0; %Branch: Name: L1_4, Typ: Line (ID: 4117b329-d609-41b6-8d18-5ed98c217161)
2	3	0.08952893	0.23380165	0.0	0	0	0	0.0	0.0	1	-360.0	360.0; %Branch: Name: L2_3, Typ: Line (ID: c39fff14-21fa-4c09-bf9f-a09b9096177c)
3	4	0.05968595	0.15586777	0.0	0	0	0	0.0	0.0	1	-360.0	360.0; %Branch: Name: L3_4, Typ: Line (ID: ac5b8882-89ea-4554-9751-8f1f40236601)
];