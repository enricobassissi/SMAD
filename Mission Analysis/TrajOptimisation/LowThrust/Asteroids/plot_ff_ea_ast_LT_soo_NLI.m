function [output, r_EA, r1,v_EA,v1] = plot_ff_ea_ast_LT_soo_NLI(x,sim,data)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDP1 = MJD01 + TOF1;

asteroid_1 = data.asteroid_names(x(5));
% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);

r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% passage at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s

r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% N REV
N_rev= x(4);

% launcher_stuff
v_inf_magn = x(6);
alpha = x(7);
beta = x(8);
v_launcher = v_inf_magn*[cos(alpha)*cos(beta); sin(alpha)*cos(beta); sin(beta)] ;
v_dep = v_EA + v_launcher;  %since parabolic escape (vinf = 0)
varr = v1; %since we want to rendez-vous

[output] = NL_interpolator( r_EA , r1 , v_dep , varr , N_rev , TOF1 , sim.M ,sim.PS.Isp ,sim );


end

