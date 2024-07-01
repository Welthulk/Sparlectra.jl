function mpc = case2
% Power flow data for 2 bus, Example 4.1 from ENTSOE_CGMES_v2.4_17Nov2014_ExplicitLoadFlowCalculation.docx
%

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 400;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	0	1.045	0	380	1	1.2	0.8;
	2	1	332 -78		0	0	0	1	0	380	1	1.2	0.8;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	400	0	1.045	400	1	400	400	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.011622160664819943	0.035263157894736843	0.0	550	550	550	0	0	1	-360	360;
	1	2	0.011622160664819943	0.035263157894736843	0.0	300	300	300 1.0021757234957092	1.77605	1	-360	360;
 ];