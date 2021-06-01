function [output, r_encounter, v_encounter,sol] = plot_ff_chesucc1_1ast_to_another_ast(x,sim,data,sol)
% setting the input times
MJD0 = sol.time_arr_ast_ADIM;
CT1 = x(1);
MJDP1 = MJD0+CT1; % dep from ast
TOF = x(2);
MJDP2 = MJDP1 + TOF;

% N REV
N_rev = x(3);

% chosing which asteroid to visit
IDP = x(4); %index of permutation, the column of the Permutation Matrix of the asteroids
% asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_1 = sol.asteroid_1;

asteroid_2 = sim.asteroid_to_fish(IDP);

% Computing position and velocity of the planets in that days
% % Departure from Earth
% MJD01_dim = MJD01*sim.TU/(3600*24);
% [kep_EA,ksun] = uplanet(MJD01_dim, 3);
% [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
% r_EA = r_EA/sim.DU;
% v_EA = v_EA/sim.DU*sim.TU;

% passage at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,sim.mu_dim); % km, km/s
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% passage at 2nd ast
MJDP2_dim = MJDP2*sim.TU/(3600*24);
[kep_ast_2] = uNEO2(MJDP2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[r2, v2] = sv_from_coe(kep_ast_2,sim.mu_dim); % km, km/s
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% v_launcher = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)] ;
% v_dep = v1 + v_launcher;  %since parabolic escape (vinf = 0)
% varr = v2; %since we want to rendezvous

[output] = NL_interpolator_of( r1 , r2 , v1 , v2 , N_rev , TOF , sim.M ,sim.PS.Isp , sim );
tol_TOF = 0.5; % meaning approx 30 days
penalty_T = 0; penalty_TOF = 0; penalty_MF_unfeasible = 0;
if max(abs(output.T_magn)) > sim.max_Available_Thrust
    penalty_T = max(abs(output.T_magn)) - sim.max_Available_Thrust;
end
if abs(output.t(end) - TOF) > tol_TOF
    penalty_TOF = abs(output.t(end) - TOF);
end

mass_fract = (output.m(1) - output.m(end))/output.m(1);
sol.mass_fract = mass_fract;

if mass_fract < 0 || mass_fract > 1
    penalty_MF_unfeasible = 100;
end

sol.obj_fun = mass_fract + penalty_MF_unfeasible + 10*penalty_T + penalty_TOF;

r_encounter.ast1 = r1;
r_encounter.ast2 = r2;

v_encounter.ast1 = v1;
v_encounter.ast2 = v2;
end

