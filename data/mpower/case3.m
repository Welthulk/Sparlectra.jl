function mpc = case3
% The network consists of 3 buses, 3 branches, and 1 generator.
%
% Network diagram:
%
%   SG->1 ---- 2 ---- 3
%
% Generator at Bus 1.
%

mpc.version = '2';
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
 1 1 0 0 0 1 1 0.8 0 0.8 0 1.2 0.8;
 2 1 0 0 0 1 1 1 0 1 0 1.2 0.8;
 3 3 0 0 0 1 1 1 0 1 0 1.2 0.8;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
 1  1.0 0.5 0.5 -0.3  0.8 100 1 0	0	0	0	0	0	0	0	0	0	0	0	0;
 3  0.0 0.0 1.5 -1.3  1.0 100 1 0	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
 1 2 0.1 0.2 0.03 0 0 0 0 0 1 -360 360;
 2 3 0.2 0.3 0.04 0 0 0 0 0 1 -360 360;
 3 1 0.15 0.25 0.02 0 0 0 0 0 0 -360 360;
];
