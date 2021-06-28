function [output, r_encounter, v_encounter, sol] = plot_ff_2RL_all_indietro2(x,sim,data, sol)
%{ 
input variable vector
x = [(1) MJD0,
     (2) TOF1,
     (3) TOF2,
     (4) TOFa,
     (5) TOFb,
     (6) NREV,
     (7) NREV2,
     (8) NREVa,
     (9) NREVb,
     (10) IDP,
     (11) IDP2_0_100,
     (12) v_inf_magn,
     (13) azimuth,
     (14) elevation,
     (15) CT1,
     (16) CT2]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDA1 = MJD01 + TOF1; % mjd2000 passage of 1st sc on ast 1
CT1 = x(15);
MJDD1 = MJDA1 + CT1;
TOF2 = x(3); % tof sc1 to 2nd asteroid
MJDA2 = MJDD1 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = x(4); % tof sc2 to 1st asteroid
MJDAa = MJD01 + TOFa; 
CTa = x(16);
MJDDa = MJDAa + CTa;
TOFb =  x(5); % tof sc2 to 2nd asteroid
MJDAb = MJDDa + TOFb; 

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if max(TOF1+CT1+TOF2,TOFa+CTa+TOFb) > max_duration
    penalty_MAX_DURATION = max(TOF1+CT1+TOF2,TOFa+CTa+TOFb) - max_duration; % 12 years max mission time 
end

% N REV1
N_rev1 = x(6);
% N REV2
N_rev2 = x(7);
% N REVa
N_reva = x(8);
% N REVb
N_revb = x(9);

% C3 launcher
v_inf_magn = x(12);
az = x(13);
elev = x(14);

%% choosing which asteroid to visit
% 1ST SPACECRAFT ASTEROID OBJECTIVES
IDP1 = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% 2ND SPACECRAFT ASTEROID OBJECTIVES
IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
asteroid_a = PermutationMatrix_2SC(IDP2,1);
asteroid_b = PermutationMatrix_2SC(IDP2,2);
% asteroid_1 = "2021DR";
% asteroid_2 = "2021JE1";
% asteroid_a = "2020VV";
% asteroid_b = "2013GM3";

%% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% ARRIVAL at 1st ast
MJDA1_dim = MJDA1*sim.TU/(3600*24);
[kep_ast_A1] = uNEO3(MJDA1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
rA1 = rA1/sim.DU;
vA1 = vA1/sim.DU*sim.TU;
% DEPARTURE at 1st ast
MJDD1_dim = MJDD1*sim.TU/(3600*24);
[kep_ast_D1] = uNEO3(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
rD1 = rD1/sim.DU;
vD1 = vD1/sim.DU*sim.TU;

% ARRIVAL at 2nd ast
MJDA2_dim = MJDA2*sim.TU/(3600*24);
[kep_ast_A2] = uNEO3(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;

% ARRIVAL at a_th ast
MJDAa_dim = MJDAa*sim.TU/(3600*24);
[kep_ast_Aa] = uNEO3(MJDAa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rAa, vAa] = sv_from_coe(kep_ast_Aa,ksun); % km, km/s
rAa = rAa/sim.DU;
vAa = vAa/sim.DU*sim.TU;
% DEPARTURE at a_th ast
MJDDa_dim = MJDDa*sim.TU/(3600*24);
[kep_ast_Da] = uNEO3(MJDDa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rDa, vDa] = sv_from_coe(kep_ast_Da,ksun); % km, km/s
rDa = rDa/sim.DU;
vDa = vDa/sim.DU*sim.TU;

% passage at b_th ast
MJDAb_dim = MJDAb*sim.TU/(3600*24);
[kep_ast_Ab] = uNEO3(MJDAb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;

%% Launcher departure variable

v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

%% NLI
Isp = sim.PS.Isp*sim.TU;
g0 = sim.g0/sim.TU^2*(1000*sim.DU);

tol_TOF = 1; % 1 TU means approx 60 days
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;

% ------------ SC1 ------------- %
%   2nd leg - Ast1 -> Ast2
[output_2] = NL_interpolator_of( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , sim.M1_end ,sim.PS.Isp , sim );
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = max(abs(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

% 1st leg - Earth -> Ast 1
M_start_1st_leg = output_2.m(1) + sim.M_pods;
[output_1] = NL_interpolator_of( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , M_start_1st_leg ,sim.PS.Isp , sim );
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg1 = max(abs(output_1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end

% ------ SC 2 ------ %
% bth leg - Ast a -> Ast b
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , sim.M2_end ,sim.PS.Isp , sim );
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = 10*max(abs(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

% a_th leg - Earth -> Ast_a
M_start_a_th_leg = output_b.m(1) + sim.M_pods;
[output_a] = NL_interpolator_of( r_EA , rAa , v_dep , vAa , N_reva , TOFa , M_start_a_th_leg ,sim.PS.Isp , sim );
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
    penalty_T_lega = 10*max(abs(output_a.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end

%% mass fractions
sol.mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

sol.mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)

sol.tot_mass_depleted_SC1 = sol.mass_depleted_Leg1+sol.mass_depleted_Leg2;
sol.mass_dry_and_pods_SC1 = output_1.m(1) - sol.tot_mass_depleted_SC1;
sol.mass_fract_SC1 = sol.tot_mass_depleted_SC1/output_1.m(1);

sol.mass_depleted_Lega = output_a.m(1) - output_a.m(end);
sol.dV_associated_Lega = -g0*Isp*log(output_a.m(end)/output_a.m(1)); % -ve*ln(m_final/m_initial)

sol.mass_depleted_Legb = output_b.m(1) - output_b.m(end);
sol.dV_associated_Legb = -g0*Isp*log(output_b.m(end)/output_b.m(1)); % -ve*ln(m_final/m_initial)

sol.tot_mass_depleted_SC2 = sol.mass_depleted_Lega+sol.mass_depleted_Legb;
sol.mass_dry_and_pods_SC2 = output_a.m(1) - sol.tot_mass_depleted_SC2;
sol.mass_fract_SC2 = sol.tot_mass_depleted_SC2/output_a.m(1);

%% obj fun
Penalty_about_feasibility_of_mass_fraction = 0;
if sol.mass_fract_SC1 < 0 || sol.mass_fract_SC1 > 1 || sol.mass_fract_SC2 < 0 || sol.mass_fract_SC2 > 1 
    Penalty_about_feasibility_of_mass_fraction = 100;
end

avg_mass_fraction = (sol.mass_fract_SC1+ sol.mass_fract_SC2)/2;
MF = max(sol.mass_fract_SC1,sol.mass_fract_SC2) + 10*abs(sol.mass_fract_SC1 - avg_mass_fraction) + ...
        10*abs(sol.mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali

sol.obj_fun = MF + penalty_MAX_DURATION + Penalty_about_feasibility_of_mass_fraction + ...
    100*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb;


%% Output encounter states
% --- new positions with actual TOFs
% ----- Computing position and velocity of the planets in that days
% 1st spacecraft characteristic times
TOF1 = output_1.t(end); % tof sc1 to 1st asteroid
MJDA1 = MJD01 + TOF1; % mjd2000 passage of 1st sc on ast 1
CT1 = x(15);
MJDD1 = MJDA1 + CT1;
TOF2 = output_2.t(end); % tof sc1 to 2nd asteroid
MJDA2 = MJDD1 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = output_a.t(end); % tof sc2 to 1st asteroid
MJDAa = MJD01 + TOFa; 
CTa = x(16);
MJDDa = MJDAa + CTa;
TOFb =  output_b.t(end); % tof sc2 to 2nd asteroid
MJDAb = MJDDa + TOFb; 

% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% ARRIVAL at 1st ast
MJDA1_dim = MJDA1*sim.TU/(3600*24);
[kep_ast_A1] = uNEO3(MJDA1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
rA1 = rA1/sim.DU;
vA1 = vA1/sim.DU*sim.TU;
% DEPARTURE at 1st ast
MJDD1_dim = MJDD1*sim.TU/(3600*24);
[kep_ast_D1] = uNEO3(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
rD1 = rD1/sim.DU;
vD1 = vD1/sim.DU*sim.TU;

% ARRIVAL at 2nd ast
MJDA2_dim = MJDA2*sim.TU/(3600*24);
[kep_ast_A2] = uNEO3(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;

% ARRIVAL at a_th ast
MJDAa_dim = MJDAa*sim.TU/(3600*24);
[kep_ast_Aa] = uNEO3(MJDAa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rAa, vAa] = sv_from_coe(kep_ast_Aa,ksun); % km, km/s
rAa = rAa/sim.DU;
vAa = vAa/sim.DU*sim.TU;
% DEPARTURE at a_th ast
MJDDa_dim = MJDDa*sim.TU/(3600*24);
[kep_ast_Da] = uNEO3(MJDDa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rDa, vDa] = sv_from_coe(kep_ast_Da,ksun); % km, km/s
rDa = rDa/sim.DU;
vDa = vDa/sim.DU*sim.TU;

% passage at b_th ast
MJDAb_dim = MJDAb*sim.TU/(3600*24);
[kep_ast_Ab] = uNEO3(MJDAb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;

r_encounter.EA = r_EA;
r_encounter.astA1 = rA1;
r_encounter.astD1 = rD1;
r_encounter.astA2 = rA2;
r_encounter.astAa = rAa;
r_encounter.astDa = rDa;
r_encounter.astAb = rAb;

v_encounter.EA = v_EA;
v_encounter.astA1 = vA1;
v_encounter.astD1 = vD1;
v_encounter.astA2 = vA2;
v_encounter.astAa = vAa;
v_encounter.astDa = vDa;
v_encounter.astAb = vAb;

%% Porcherie
t_span_CT1 = linspace(output_1.t(end),output_1.t(end)+CT1,sim.n_sol);
t_span_CTa = linspace(output_a.t(end),output_a.t(end)+CTa,sim.n_sol);
% t_span_CT1 = linspace(0,CT1,sim.n_sol);
% t_span_CTa = linspace(0,CTa,sim.n_sol);
mCT1 = ones(sim.n_sol,1).*output_1.m(end);
mCTa = ones(sim.n_sol,1).*output_a.m(end);
TCT1 = zeros(sim.n_sol,3); % .*output_a.Thrust(end)
TCTa = zeros(sim.n_sol,3);

%% Output
output.t_SC1            = [output_1.t; t_span_CT1'; t_span_CT1(end)+output_2.t];
output.t_SC2            = [output_a.t; t_span_CTa'; t_span_CTa(end)+output_b.t];
output.m_SC1            = [output_1.m; mCT1; output_2.m];
output.m_SC2            = [output_a.m; mCTa; output_b.m];
output.Thrust_SC1       = [output_1.Thrust; TCT1; output_2.Thrust];
output.T_magn_SC1       = sqrt(output.Thrust_SC1(:,1).^2 + output.Thrust_SC1(:,3).^2);
output.Thrust_SC2       = [output_a.Thrust; TCTa; output_b.Thrust];
output.T_magn_SC2       = sqrt(output.Thrust_SC2(:,1).^2 + output.Thrust_SC2(:,3).^2);
output.a_SC1            = [output_1.a; output_2.a];
output.a_SC2            = [output_a.a; output_b.a];
output.r.leg1       = output_1.r; 
output.r.leg2       = output_2.r;
output.r.lega       = output_a.r;
output.r.legb       = output_b.r;
output.theta.leg1   = output_1.theta; 
output.theta.leg2   = output_2.theta;
output.theta.lega   = output_a.theta;
output.theta.legb   = output_b.theta;
output.z.leg1       = output_1.z;
output.z.leg2       = output_2.z;
output.z.lega       = output_a.z;
output.z.legb       = output_b.z;
output.Href.leg1    = output_1.Href;
output.Href.leg2    = output_2.Href;
output.Href.lega    = output_a.Href;
output.Href.legb    = output_b.Href;

output.t1           = output_1.t;
output.CT1          = t_span_CT1';
output.t2           = output_2.t;
output.ta           = output_a.t;
output.CTa          = t_span_CTa';
output.tb           = output_b.t;

%% --
sol.T_1 = output_1.Thrust;
sol.T_2 = output_2.Thrust;
sol.T_a = output_a.Thrust;
sol.T_b = output_b.Thrust;

end

