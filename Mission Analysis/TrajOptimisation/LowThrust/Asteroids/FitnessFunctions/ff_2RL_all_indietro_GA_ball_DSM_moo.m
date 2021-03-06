function obj_fun = ff_2RL_all_indietro_GA_ball_DSM_moo(x,sim,data)
%{ 

input variable vector
% x = [Departure dates (1)
% TOF_DSM1 (2), from earth dep to DSM time
% rD (3)
% thetaD (4)
% TOFGA (5), from DSM1 to GA
% azimuth out (POST-GRAVITY ASSIST)(6)
% elevation (POST-GRAVITY ASSIST)(7)
% TOF1 (8)
% N REV 1 (9)
% Coasting time 1 (10)
% TOF2 (11)
% N REV 2 (12)
% TOF_DSM2 (13) from earth dep to DSM time -> determined by thetaD2
% rD2 (14)
% thetaD2 (15)
% TOFGAa (16)
% azimuth out (POST-GRAVITY ASSIST)(17)
% elevation (POST-GRAVITY ASSIST)(18)
% TOFa (19)
% N REV a (20)
% Coasting time a (21)
% TOFb (22)
% N REV b (23)
% -- ID Permutation (24)
% ID Permutation 2 (25)]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc
% 1st spacecraft characteristic times
TOF_DSM1 = x(2); % It's the time from Earth to DSM
TOFGA_1 = x(5); % DSM to GA
MJDPGA_1 = MJD01 + TOF_DSM1 + TOFGA_1; % passage time GA
TOF1 = x(8); % tof sc1 GA to Ast1 
MJDA1 = MJDPGA_1 + TOF1; % mjd2000 arrive of 1st sc on ast 1
CT1 = x(10);
MJDD1 = MJDA1 + CT1; % departure from 1st asteroid
TOF2 = x(11); % tof sc1 from ast 1 to 2nd asteroid
MJDA2 = MJDD1 + TOF2; % mjd2000 passage of SC1 on ast 2

% ---- SC2 times
TOF_DSM2 = x(13);
TOFGA_a = x(16); % DSM to GA
MJDPGA_a = MJD01 + TOF_DSM2 + TOFGA_a; % passage time GA
TOFa = x(19); % tof sc1 GA to Ast1 
MJDAa = MJDPGA_a + TOFa; % mjd2000 arrive of 1st sc on ast 1
CTa = x(21);
MJDDa = MJDAa + CTa; % departure from 1st asteroid
TOFb = x(22); % tof sc1 from ast 1 to 2nd asteroid
MJDAb = MJDDa + TOFb; % mjd2000 passage of SC1 on ast 2

% -- N REV1
N_rev1 = x(9);
% N REV2
N_rev2 = x(12);
% N REVa
N_reva = x(20);
% N REVb
N_revb = x(23);

% -- GA stuff SC1
% OUT
az_out_GA1 = x(6);
el_out_GA1 = x(7);

% GA stuff SC2
% OUT
az_out_GAa = x(17);
el_out_GAa = x(18);

% -- DSM SC1
rD_mag1 = x(3);
thetaD1 = x(4);

% -- DSM SC2
rD_mag2 = x(14);
thetaD2 = x(15);

%% choosing which asteroid to visit
% 1ST SPACECRAFT ASTEROID OBJECTIVES
IDP1 = x(24); %index of permutation, the column of the Permutation Matrix of the asteroids
% IDP1 = ceil(x(28)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% % 2ND SPACECRAFT ASTEROID OBJECTIVES
IDP_temp_2 = x(25); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
% IDP_temp_2 = ceil(x(29)); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
asteroid_a = PermutationMatrix_2SC(IDP2,1);
asteroid_b = PermutationMatrix_2SC(IDP2,2);

%% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% -- DSM
iEA = kep_EA(3); % same plane as earth at beginning
rDSM1 = rD_mag1*[cos(thetaD1)*cos(iEA); sin(thetaD1)*cos(iEA); sin(iEA)];

% GA SC1
MJDPGA1_dim = MJDPGA_1*sim.TU/(3600*24);
[kep_GA1,ksun] = uplanet(MJDPGA1_dim, sim.ID_FLYBY);
[r_GA1, v_GA1] = sv_from_coe(kep_GA1,ksun);
r_GA1 = r_GA1/sim.DU; 
v_GA1 = v_GA1/sim.DU*sim.TU;

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

% --- SC positions
% -- DSM
iEA = kep_EA(3); % same plane as earth at beginning
rDSM2 = rD_mag2*[cos(thetaD2)*cos(iEA); sin(thetaD2)*cos(iEA); sin(iEA)];

% GA SC2
MJDPGAa_dim = MJDPGA_a*sim.TU/(3600*24);
[kep_GAa,ksun] = uplanet(MJDPGAa_dim, 3);
[r_GAa, v_GAa] = sv_from_coe(kep_GAa,ksun);
r_GAa = r_GAa/sim.DU;
v_GAa = v_GAa/sim.DU*sim.TU;

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

%% FUEL CONSUMPTIONS %%
tol_TOF = 1; % 1 TU means approx 60 days
penalty_dv_GA1 = 0; penalty_dv_GAa = 0;
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;

Isp_secondary_prop = sim.PS.Isp_secondary_prop*sim.TU; % ADN
g0 = 9.81;

%% Both SC togheter
% ballistic stuff before GA
% DV calculation with lambert
% Earth -> DSM
[~,~,~,~,V_EA_DSM_dep,V_EA_DSM_arr,~,~] = lambertMR(r_EA,rDSM1,TOF_DSM1,sim.mu,0,0,0,0);
V_EA_DSM_dep = V_EA_DSM_dep';
V_EA_DSM_arr = V_EA_DSM_arr';
sqrt_C3_requested = sqrt((V_EA_DSM_dep(1)-v_EA(1))^2+(V_EA_DSM_dep(2)-v_EA(2))^2+(V_EA_DSM_dep(3)-v_EA(3))^2);
% C3 vettore ?? solo V_EA_DSM_dep - v_EA
if sqrt_C3_requested < sqrt(sim.C3_max) % vinf that the launcher can give max 
    penalty_dV_C3 = 0;
else
    penalty_dV_C3 = 100*(sqrt_C3_requested - sqrt(sim.C3_max)); % penalty like, but not discard a priori
end

%% SC 1
% 2nd leg - Ast1 -> Ast2
[output_2] = NL_interpolator_of( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , sim.M2_end ,sim.PS.Isp , sim );
%     if ~isnan(output_2.Thrust(1,1)) && abs(output_2.t(end) - TOF2) < tol_TOF  % if is not nan -> it's a valid solution
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = max(abs(output_2.T_magn));
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

% DSM -> GA
[~,~,~,~,V_DSM_GA_dep,V_DSM_GA_arr,~,~] = lambertMR(rDSM1,r_GA1,TOFGA_1,sim.mu,0,0,0,0);
V_DSM_GA_dep = V_DSM_GA_dep';
V_DSM_GA_arr = V_DSM_GA_arr';
dV_DSM_GA = sqrt((V_DSM_GA_dep(1)-V_EA_DSM_arr(1))^2+(V_DSM_GA_dep(2)-V_EA_DSM_arr(2))^2+(V_DSM_GA_dep(3)-V_EA_DSM_arr(3))^2);
if dV_DSM_GA < sim.dV_DSM_max % vinf that the launcher can give max 
    pen_dV_DSM = 0;
else
    pen_dV_DSM = 100*abs(dV_DSM_GA - sim.dV_DSM_max); % penalty like, but not discard a priori
end

% 1st leg - GA -> Ast 1
M_end_1st_leg = output_2.m(1) + sim.M_pods; % yes pods here, released on Ast1
v_inf_minus = V_DSM_GA_arr - v_GA1;
v_inf_plus = norm(v_inf_minus)*[cos(el_out_GA1)*cos(az_out_GA1); cos(el_out_GA1)*sin(az_out_GA1); sin(el_out_GA1)];
v_dep_GA1 = v_inf_plus + v_GA1; 
[output_1] = NL_interpolator_of( r_GA1 , rA1 , v_dep_GA1 , vA1 , N_rev1 , TOF1 , M_end_1st_leg ,sim.PS.Isp , sim );
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg1 = max(abs(output_1.T_magn));
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end

% Gravity Assist
if sim.ID_FLYBY == 3
    RPlanet_flyby = astroConstants(23); % Radius_Earth, km
    muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
    R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
    R_SOI_PL = 0.929*1e6; % km
elseif sim.ID_FLYBY == 4
    RPlanet_flyby = astroConstants(24); % Radius_mars, km
    muPlanet_flyby = astroConstants(14); % mu mars, km^3/s^2
    R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
    R_SOI_PL = 0.578*1e6; %km
end
v_arr_GA1_dim = V_DSM_GA_arr.*sim.DU./sim.TU;
v_dep_GA1_dim = v_dep_GA1.*sim.DU./sim.TU;
[delta_v_p_1,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGA1_dim, v_arr_GA1_dim, v_dep_GA1_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_1), 'Not found')
    penalty_dv_GA1 = 10;
end

Mass_at_start_SC1 = output_1.m(1)*exp(dV_DSM_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));

%% SC 2
% Earth -> DSM2
[~,~,~,~,V_EA_DSM2_dep,V_EA_DSM2_arr,~,~] = lambertMR(r_EA,rDSM2,TOF_DSM2,sim.mu,0,0,0,0);
V_EA_DSM2_dep = V_EA_DSM2_dep';
V_EA_DSM2_arr = V_EA_DSM2_arr';
% i will pay the difference between actual c3 applied, the one coming out
% from the optimisation of SC1, and the one required for this leg
dV_diff_wrt_C3_SC1 = sqrt((V_EA_DSM2_dep(1)-V_EA_DSM_dep(1))^2+(V_EA_DSM2_dep(2)-V_EA_DSM_dep(2))^2+(V_EA_DSM2_dep(3)-V_EA_DSM_dep(3))^2);

% bth leg - Asta -> Astb
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , sim.M2_end ,sim.PS.Isp , sim );
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = max(abs(output_b.T_magn));
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

% DSM -> GA
[~,~,~,~,V_DSM2_GA_dep,V_DSM2_GA_arr,~,~] = lambertMR(rDSM2,r_GAa,TOFGA_a,sim.mu,0,0,0,0);
V_DSM2_GA_dep = V_DSM2_GA_dep';
V_DSM2_GA_arr = V_DSM2_GA_arr';
dV_DSM2_GA = sqrt((V_DSM2_GA_dep(1)-V_EA_DSM2_arr(1))^2+(V_DSM2_GA_dep(2)-V_EA_DSM2_arr(2))^2+(V_DSM2_GA_dep(3)-V_EA_DSM2_arr(3))^2);
if dV_DSM2_GA < sim.dV_DSM_max % vinf that the launcher can give max 
    pen_dV_DSM2 = 0;
else
    pen_dV_DSM2 = 100*abs(dV_DSM2_GA - sim.dV_DSM_max); % penalty like, but not discard a priori
end

% ath leg - GA -> Ast a
M_end_bth_leg = output_b.m(1) + sim.M_pods; % yes pods here, released on Ast1
v_inf_minus_a = V_DSM2_GA_arr - v_GAa;
v_inf_plus_a = norm(v_inf_minus_a)*[cos(el_out_GAa)*cos(az_out_GAa); cos(el_out_GAa)*sin(az_out_GAa); sin(el_out_GAa)];
v_dep_GAa = v_inf_plus_a + v_GAa; 
[output_a] = NL_interpolator_of( r_GAa , rAa , v_dep_GAa , vAa , N_reva , TOFa , M_end_bth_leg ,sim.PS.Isp , sim );
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
    penalty_T_lega = max(abs(output_a.T_magn));
end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end

% Gravity Assist
if sim.ID_FLYBY == 3
    RPlanet_flyby = astroConstants(23); % Radius_Earth, km
    muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
    R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
    R_SOI_PL = 0.929*1e6; % km
elseif sim.ID_FLYBY == 4
    RPlanet_flyby = astroConstants(24); % Radius_mars, km
    muPlanet_flyby = astroConstants(14); % mu mars, km^3/s^2
    R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
    R_SOI_PL = 0.578*1e6; %km
end
v_arr_GAa_dim = V_DSM2_GA_arr.*sim.DU./sim.TU;
v_dep_GAa_dim = v_dep_GAa.*sim.DU./sim.TU;
[delta_v_p_a,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGAa_dim, v_arr_GAa_dim, v_dep_GAa_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_a), 'Not found')
    penalty_dv_GAa = 10;
end

Mass_at_start_SC2 = output_a.m(1)*exp(dV_DSM2_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));

%% masses and obj function
% SC1
mass_depleted_DSM1 = Mass_at_start_SC1 - output_1.m(1);
mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);

tot_mass_depleted_SC1 = mass_depleted_DSM1+mass_depleted_Leg1+mass_depleted_Leg2;
mass_fract_SC1 = tot_mass_depleted_SC1/Mass_at_start_SC1;

% SC2
mass_depleted_DSM2 = Mass_at_start_SC2 - output_a.m(1);
mass_depleted_Lega = output_a.m(1) - output_a.m(end);
mass_depleted_Legb = output_b.m(1) - output_b.m(end);

tot_mass_depleted_SC2 = mass_depleted_DSM2+mass_depleted_Lega+mass_depleted_Legb;
mass_fract_SC2 = tot_mass_depleted_SC2/Mass_at_start_SC2;

penalty_MF_unfeasible = 0;
if mass_fract_SC1 <= 0 || mass_fract_SC1 >= 1 || mass_fract_SC2 <= 0 || mass_fract_SC2 >= 1
    penalty_MF_unfeasible = 100;
end

avg_mass_fraction = (mass_fract_SC1+mass_fract_SC2)/2;
MF = max(mass_fract_SC1,mass_fract_SC2) + abs(mass_fract_SC1 - avg_mass_fraction) + ...
        abs(mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if max(TOF_DSM1+TOFGA_1+TOF1+CT1+TOF2,TOF_DSM2+TOFGA_a+TOFa+CTa+TOFb) > max_duration
    penalty_MAX_DURATION = max(TOF_DSM1+TOFGA_1+TOF1+CT1+TOF2,TOF_DSM2+TOFGA_a+TOFa+CTa+TOFb) - max_duration; % 12 years max mission time 
end

% 2 obj function
% 1: related to mass fraction and dV + all unfeasibility stuff
% 2: related to T and excess of thrust + all unfeasibility stuff

% obj_fun(1) = 20*(dV_DSM_GA + dV_DSM2_GA) + penalty_dV_C3 + 100*dV_diff_wrt_C3_SC1 + ...
%     pen_dV_DSM + pen_dV_DSM2 + penalty_dv_GA1 + penalty_dv_GAa + ...
%     MF + penalty_MAX_DURATION + penalty_MF_unfeasible + ...
%     penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb;
obj_fun(1) = 20*(dV_DSM_GA + dV_DSM2_GA) + penalty_dV_C3 + 100*dV_diff_wrt_C3_SC1 + ...
    pen_dV_DSM + pen_dV_DSM2 + penalty_dv_GA1 + penalty_dv_GAa + ...
    MF + penalty_MAX_DURATION + penalty_MF_unfeasible;

% obj_fun(2) = 20*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
%     penalty_MAX_DURATION + penalty_MF_unfeasible + ...
%     penalty_dv_GA1 + penalty_dv_GAa + ...
%     penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb;
obj_fun(2) = 30*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
    penalty_MAX_DURATION + penalty_MF_unfeasible + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb;

end