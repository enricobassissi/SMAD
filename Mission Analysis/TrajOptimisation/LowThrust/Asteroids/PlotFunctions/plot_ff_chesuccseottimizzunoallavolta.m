function [output, r_encounter, v_encounter,sol] = plot_ff_chesuccseottimizzunoallavolta(x,sim,data,sol)
% setting the input times
MJD01 = x(1);
TOF1 = x(2);
MJDP1 = MJD01 + TOF1;

% N REV
N_rev = x(3);

% C3 launcher
v_inf_magn = x(4);
az = x(5);
el = x(6);

% chosing which asteroid to visit
IDP = x(7); %index of permutation, the column of the Permutation Matrix of the asteroids
% asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_1 = sim.asteroid_to_fish(IDP); % data.asteroid_names(IDP); 
% asteroid_2 = data.PermutationMatrix(IDP,2);
% asteroid_3 = data.PermutationMatrix(IDP,3);
% asteroid_4 = data.PermutationMatrix(IDP,4);

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

v_launcher = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)] ;
v_dep = v_EA + v_launcher;  %since parabolic escape (vinf = 0)
% varr = v2; %since we want to rendezvous

[output] = NL_interpolator_of( r_EA , r1 , v_dep , v1 , N_rev , TOF1 , sim.M ,sim.PS.Isp , sim );
tol_TOF = 0.5; % meaning approx 30 days
penalty_T = 0; penalty_TOF = 0; penalty_MF_unfeasible = 0;
if max(abs(output.T_magn)) > sim.max_Available_Thrust
    penalty_T = max(abs(output.T_magn)) - sim.max_Available_Thrust;
end
if abs(output.t(end) - TOF1) > tol_TOF
    penalty_TOF = abs(output.t(end) - TOF1);
end

mass_fract = (output.m(1) - output.m(end))/output.m(1);
sol.mass_fract = mass_fract;

if mass_fract < 0 || mass_fract > 1
    penalty_MF_unfeasible = 100;
end

sol.obj_fun = mass_fract + penalty_MF_unfeasible + 10*penalty_T + penalty_TOF;

r_encounter.EA = r_EA;
r_encounter.ast1 = r1;

v_encounter.EA = v_EA;
v_encounter.ast1 = v1;

end

