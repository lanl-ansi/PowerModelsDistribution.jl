% This case was designed for testing the voltage balance constraints

% 3 coupled radial networks
% no shunts
% This case contains a load in phase 1 at bus 2,
% a cheap generator at bus 3, and an expensive one at bus 2
% The optimal dispatch will use the generator at bus 3 as much as possible,
% being constrained by the imposed balance constraint

function tppmc = case5_i_r_a
tppmc.version = '1'

tppmc.baseMVA = 100.0;
tppmc.baseKV = 230.0;

%% bus data
%	bus_i	type	vmin_1	vmax_1	vmin_2	vmax_2	vmin_3	vmax_3	vm_1	va_1	vm_2	va_2	vm_3	va_3
tppmc.bus = [
	1	 3	 0.9999	 1.0001	 0.9999	 1.0001	 0.9999	 1.0001	 1.00000	 2.80377	 1.00000	 2.80377	 1.00000	 2.80377;
	2	 1	 0.90	 1.10	 0.90	 1.10	 0.90	 1.10	 1.08407	-0.73465	 1.08407	-0.73465	 1.08407	-0.73465;
	3	 1	 0.90	 1.10	 0.90	 1.10	 0.90	 1.10	 1.00000	-0.55972	 1.00000	-0.55972	 1.00000	-0.55972;
];

%% load data
%	load_bus	pd_1	qd_1	pd_2	qd_2	pd_3	qd_3	status
tppmc.load = [
	2	 100.0	 100.0	 0.0	 0.0	 0.0	 0.0	 1;
];

%% generator data
%	gen_bus	pmin_1	pmax_1	qmin_1	qmax_1	pmin_2	pmax_2	qmin_2	qmax_2	pmin_3	pmax_3	qmin_3	qmax_3	pg_1	qg_1	pg_2	qg_2	pg_3	qg_3	status
tppmc.gen = [
	2	 0.0	  40.0	  -30.0	  30.0	 0.0	  40.0	  -30.0	  30.0	 0.0	  40.0	  -30.0	  30.0	  40.0	  30.0	  40.0	  30.0	  40.0	  30.0	 1;
	3	 0.0	 170.0	 -127.5	 127.5	 0.0 	 170.0	 -127.5	 127.5	 0.0 	 170.0	 -127.5	 127.5	 170.0	 127.5	 170.0	 127.5	 170.0	 127.5	 1;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
tppmc.gencost = [
	2	 0.0	 0.0	 3	 0.0	 100.0	 0.0;
	2	 0.0	 0.0	 3	 0.0	 10.0	 0.0;
];

%% branch data
%	f_bus	t_bus	r_11	x_11	r_12	x_12	r_13	x_13	r_22	x_22	r_23	x_23	r_33	x_33   b_1    b_2    b_3    rate_a	rate_b	rate_c	angmin	angmax	status
tppmc.branch = [
	1	 2	 0.00281	 0.0281	 0.0001   0.0002	 0.0001   0.0002	 0.00281	 0.0281	 0.0001   0.0002	 0.00281	 0.0281    0.0    0.0    0.0    400	 400	 400	 -30.0	 30.0	 1;
	2	 3	 0.0108	 0.108	 0.001   0.002	 0.001   0.002	 0.0108	 0.108	 0.001   0.002	 0.0108	 0.108    0.0    0.0    0.0    426	 426	 426	 -30.0	 30.0	 1;
];

%% bus names
tppmc.bus_name = {
	'Bus A -';
	'Bus B @';
	'Bus C $';
};
