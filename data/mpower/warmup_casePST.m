function mpc = warmup_casePST
% WARMUP_CASEPST  PST example network from the Power Grid Model workshop
% (Generic Branch Representation for Line, Transformer and PST), rebuilt
% as a standard MATPOWER case.
%
% Conventions (decided against the Power Grid Model reference anchor):
% - baseMVA = 100; branch impedances converted with the to-side voltage
%   base (230 kV: Z_base = 529 Ohm).
% - PST23 (B2->B3): theta = -0.1 rad from the source model. This file is
%   written for the Sparlectra import DEFAULTS (matpower_import.shift_unit
%   = deg, shift_sign = 1.0, ratio = normal); the SHIFT sign below is the
%   one that reproduces the source anchor (L34 loading 81 percent at
%   -0.1 rad, 53 percent at 0.0 rad).
% - G_L34 = 4e-6 S is the parallel CONDUCTANCE of cable L34 in the source
%   model (dielectric/insulation losses, about 0.21 MW total). Standard
%   MATPOWER branches have no parallel-g column, so DECISION: split it as
%   bus GS onto the two cable ends B3 and B4 (0.1058 MW each at 1.0 pu).
%   In power-flow results this shows up as a small ACTIVE shunt draw
%   Ps = Gs*V^2 at B3/B4 (about 0.109 and 0.105 MW) with Qs empty
%   (Bs = 0): these shunts are loss elements, not capacitor banks. The
%   cable's large CAPACITIVE charging sits in the branch itself
%   (b = 2.116 pu) and appears as about -207 MVar on the L34 branch row.
% - Branch order: 1 L27, 2 L28, 3 L45, 4 L46, 5 L57, 6 L34, 7 T12, 8 PST23.
% - B9 (KingsCross) is a busbar section linked to B2 by a closed
%   impedance-less coupler (mpc.sparlectra.links); L28 hangs off B9.
%
% Single-line diagram (base case; ratings in MVA, loads in MW):
%
%   B1 15.5 kV (slack, 1.0 pu)
%    |
%   T12 (1000)
%    |
%   B2 230 kV ---------L27 (500)--------- B7 230 kV [200 MW]
%    |  \                                  |
%    |   --link-- B9 230 kV (busbar)      |
%    |             \                       |
%    |              L28 (500)--- B8 230 kV |  [200 MW at B8]
%    |                                     |
%   PST23 (1000, theta = -0.1 rad)     L57 (250)
%    |                                     |
%   B3 230 kV                              |
%    |                                     |
%   L34 cable (500, b = 2.116 pu)          |
%    |                                     |
%   B4 230 kV ------L45 (500)--------- B5 230 kV [100 MW]
%       \
%        ----------L46 (500)--------- B6 230 kV [100 MW]
%
% The busbar split at B2 is now part of this file: B9 (KingsCross) is a
% busbar section joined to B2 by the closed impedance-less coupler in
% mpc.sparlectra.links, and L28 feeds B8 from B9. The remaining expansion
% stage (case_udo_build.jl, not part of this file) hangs a 110-kV mesh
% (B10/B11) with two OLTC transformers under B3/B4 plus a switched MSR.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	15.5	1	1.1	0.9;
	2	1	0	0	0	0	1	1	0	230	1	1.1	0.9;
	3	1	0	0	0.1058	0	1	1	0	230	1	1.1	0.9;
	4	1	0	0	0.1058	0	1	1	0	230	1	1.1	0.9;
	5	1	100	0	0	0	1	1	0	230	1	1.1	0.9;
	6	1	100	0	0	0	1	1	0	230	1	1.1	0.9;
	7	1	200	0	0	0	1	1	0	230	1	1.1	0.9;
	8	1	200	0	0	0	1	1	0	230	1	1.1	0.9;
	9	1	0	0	0	0	1	1	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	600	0	9999	-9999	1.0	100	1	2000	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	2	7	0.00378072	0.01512287	0	500	500	500	0	0	1	-360	360;
	9	8	0.00378072	0.01512287	0	500	500	500	0	0	1	-360	360;
	4	5	0.00378072	0.01512287	0	500	500	500	0	0	1	-360	360;
	4	6	0.00378072	0.01512287	0	500	500	500	0	0	1	-360	360;
	5	7	0.00378072	0.01512287	0	250	250	250	0	0	1	-360	360;
	3	4	0.00869565	0.00642722	2.116	500	500	500	0	0	1	-360	360;
	1	2	0.00094518	0.00378072	0	1000	1000	1000	1	0	1	-360	360;
	2	3	0	0.00756144	0	1000	1000	1000	1	-5.7295780	1	-360	360;
];

%% bus names (imported with matpower_import.apply_bus_names)
%% Harry Potter naming: B1 Hogwarts, B2 Hogsmeade, B3 DiagonAlley,
%% B4 Gringotts, B5 Azkaban, B6 GodricsHollow, B7 TheBurrow, B8 PrivetDrive,
%% B9 KingsCross (busbar section of the Hogsmeade station)
mpc.bus_name = {
	'Hogwarts';
	'Hogsmeade';
	'DiagonAlley';
	'Gringotts';
	'Azkaban';
	'GodricsHollow';
	'TheBurrow';
	'PrivetDrive';
	'KingsCross';
};

%% branch names (imported with matpower_import.apply_branch_names)
mpc.branch_name = {
	'L27';
	'L28';
	'L45';
	'L46';
	'L57';
	'L34';
	'T12';
	'PST23';
};

%% Sparlectra extension: impedance-less busbar coupler (BusLink).
%% B9 KingsCross is a busbar section of the Hogsmeade station; the closed
%% link merges it galvanically with B2, so L28 (B9->B8) still carries its
%% 200 MW and the power-flow solution is unchanged. Opening the link is
%% deliberate test material for topology validation.
mpc.sparlectra = struct();
mpc.sparlectra.format_version = 1;
mpc.sparlectra.links = [
	2	9	1;
];

%% Tap-changer nameplate data (columns: branch, tap_step [fraction/step],
%% tap_min_step, tap_max_step, tap_current_step, phase_step_deg,
%% phase_min_step, phase_max_step, phase_current_step, psi_deg,
%% phase_du_step). The mpc.branch TAP/SHIFT columns stay the NEUTRAL
%% position; current steps move the live position off it. T12 (branch 7):
%% machine transformer, +-2 steps of 2.5 percent, at neutral. PST23
%% (branch 8): ADDITIONAL-VOLTAGE stepper (Delta-u PST, tap_step 0 = no
%% ratio tap changer, phase_step_deg 0): each step adds du = 0.01 pu of
%% additional voltage in direction psi = 90 deg, +-10 steps around the
%% neutral shift of -0.1 rad; the shift angle FOLLOWS from the cascade
%% (atan), it is not the mechanical grid. Standing at neutral.
mpc.sparlectra.tap_changers = [
	7	0.025	-2	2	0	0	0	0	0	0	0;
	8	0	0	0	0	0	-10	10	0	90	0.01;
];
