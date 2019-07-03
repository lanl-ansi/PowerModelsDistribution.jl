%
% Tests when a large shunt needs to shed to ensure AC feasiblity
% removing line 3-2 requires shunt shedding for feasibilty
%

function mpc = case3_s_ml
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 00.0	 0.0	 0.0	 1	    1.00000	   0.00000	 240.0	 1	    1.10000	    0.90000;
	2	 2	 100.0	 50.0	 1.0	-30.0	 1	    1.00000	   0.00000	 240.0	 1	    1.10000	    0.90000;
	3	 2	 100.0	 50.0	 0.0	-30.0	 1	    1.00000	   0.00000	 240.0	 1	    1.10000	    0.90000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 0.000	 0.000	 1000.0	 -1000.0	 1.00000	 100.0	 1	 100.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	3	 0.000	 0.000	 1000.0	 -1000.0	 1.00000	 100.0	 1	 50.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   0.000000	  10.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	   1.000000	   0.000000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 3	 0.065	 0.62	 0.0	 30.0	 0.0	 0.0	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 2	 0.025	 0.75	 0.0	 30.0	 0.0	 0.0	 0.0	 0.0	 0	 -30.0	 30.0;
	1	 2	 0.042	 0.9	 0.0	 60.0	 0.0	 0.0	 0.0	 0.0	 1	 -30.0	 30.0;
];
