% This file was generated by Sparlectra exportMatPower.jl
function mpc = case5
%% MATPOWER Case Format : Version 2
mpc.version = '2';
%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
mpc.bus = [
%bus	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
1	1	1.0	2.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: Bus_1_110 (ID: #Bus_1_110#)
2	1	1.0	2.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: Bus_2_110 (ID: #Bus_2_110#)
3	1	1.0	2.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: Bus_3_110 (ID: #Bus_3_110#)
4	1	0.0	0.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: Bus_4_110 (ID: #Bus_4_110#)
5	3	0.0	0.0	0.0	0.0	1	1.0	0.0	110.0	1	1.1	0.9; %Bus: Bus_5_110* (slack) (ID: #Bus_5_110*#)
];
%% generator data
mpc.gen = [
%bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
1	1.1	2.0	100.0	-100.0	1.0	100.0	1	100.0	-100.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0; %Gen: Generator_Gen_110)
5	0.0	0.0	100.0	-100.0	1.0	100.0	1	100.0	-100.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0; %Gen: Slack-Generator (ID: 5)
];
%% branch data
mpc.branch = [
%fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
1	2	0.00165289	0.00322314	0.0	0	0	0	0.0	0.0	1	-360	360; %Branch: B_ACL_110_1_2 (ID: #B_ACL_110_1_2#1)
1	3	0.00165289	0.00322314	0.0	0	0	0	0.0	0.0	1	-360	360; %Branch: B_ACL_110_1_3 (ID: #B_ACL_110_1_3#2)
2	4	0.00165289	0.00322314	0.0	0	0	0	0.0	0.0	1	-360	360; %Branch: B_ACL_110_2_4 (ID: #B_ACL_110_2_4#3)
3	4	0.00165289	0.00322314	0.0	0	0	0	0.0	0.0	1	-360	360; %Branch: B_ACL_110_3_4 (ID: #B_ACL_110_3_4#4)
4	5	0.00165289	0.00322314	0.0	0	0	0	0.0	0.0	1	-360	360; %Branch: B_ACL_110_4_5 (ID: #B_ACL_110_4_5#5)
];
