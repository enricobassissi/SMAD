function [output, r_encounter, v_encounter, sol] = plot_ff_2RL_GA_unp_anna2(x,sim,data,sol,power_propulsion_data)
% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF_DSM1 = x(2); % It's the time from Earth to DSM
TOFGA_1 = x(5); % DSM to GA
MJDPGA_1 = MJD01 + TOF_DSM1 + TOFGA_1; % passage time GA
TOF1 = x(11); % tof sc1 GA to Ast1 
MJDA1 = MJDPGA_1 + TOF1; % mjd2000 arrive of 1st sc on ast 1
CT1 = x(13);
MJDD1 = MJDA1 + CT1; % departure from 1st asteroid
TOF2 = x(14); % tof sc1 from ast 1 to 2nd asteroid
MJDA2 = MJDD1 + TOF2; % mjd2000 passage of SC1 on ast 2

% -- N REV1
N_rev1 = round(x(12));
% N REV2
N_rev2 = round(x(15));
% N REVa
N_reva = round(x(24));
% N REVb
N_revb = round(x(27));

% -- GA stuff SC1
% IN
v_inf_magn_GA1 = x(6);
az_in_GA1 = x(7);
el_in_GA1 = x(8);
% OUT
az_out_GA1 = x(9);
el_out_GA1 = x(10);

% GA stuff SC2
% IN
v_inf_magn_GAa = x(18);
az_in_GAa = x(19);
el_in_GAa = x(20);
% OUT
az_out_GAa = x(21);
el_out_GAa = x(22);

% -- DSM SC1
rD_mag1 = x(3);
thetaD1 = x(4);

% -- DSM SC2
% rD_mag2 = x(17);
thetaD2 = x(16);

%% choosing which asteroid to visit
% 1ST SPACECRAFT ASTEROID OBJECTIVES
IDP1 = round(x(28)); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% % 2ND SPACECRAFT ASTEROID OBJECTIVES
IDP_temp_2 = round(x(29)); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
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
[kep_ast_A1] = uNEO2(MJDA1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
rA1 = rA1/sim.DU;
vA1 = vA1/sim.DU*sim.TU;
% DEPARTURE at 1st ast
MJDD1_dim = MJDD1*sim.TU/(3600*24);
[kep_ast_D1] = uNEO2(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
rD1 = rD1/sim.DU;
vD1 = vD1/sim.DU*sim.TU;

% ARRIVAL at 2nd ast
MJDA2_dim = MJDA2*sim.TU/(3600*24);
[kep_ast_A2] = uNEO2(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;



%% Launcher departure variable
% v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
% v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

% Gravity assist SC1
v_arr_GA1 = v_GA1 + v_inf_magn_GA1*[cos(el_in_GA1)*cos(az_in_GA1); cos(el_in_GA1)*sin(az_in_GA1); sin(el_in_GA1)];

% Gravity assist SC2
% will be done later in SC2 stuff computations

%% FUEL CONSUMPTIONS %%
tol_TOF = 1; % 1 TU means approx 60 days
penalty_dv_GA1 = 0; penalty_dv_GAa = 0;
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;

Isp_secondary_prop = 255; % ADN
g0 = 9.81;

%% Both SC togheter
% ballistic stuff before GA
% DV calculation with lambert
% Earth -> DSM
[~,~,~,~,V_EA_DSM_dep,V_EA_DSM_arr,~,~] = lambertMR(r_EA,rDSM1,TOF_DSM1,sim.mu,0,0,0,0);
V_EA_DSM_dep = V_EA_DSM_dep';
V_EA_DSM_arr = V_EA_DSM_arr';
sqrt_C3_requested = norm(V_EA_DSM_dep - v_EA);
% C3 vettore è solo V_EA_DSM_dep - v_EA
if sqrt_C3_requested < sqrt(sim.C3_max) % vinf that the launcher can give max 
    penalty_dV_C3 = 0;
else
    penalty_dV_C3 = 100*(sqrt_C3_requested - sqrt(sim.C3_max)); % penalty like, but not discard a priori
end

%% SC 1
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
v_dep_GA1 = v_GA1 + v_inf_magn_GA1*[cos(el_out_GA1)*cos(az_out_GA1); cos(el_out_GA1)*sin(az_out_GA1); sin(el_out_GA1)];
Mass_at_beginning_LT = sim.M1*exp(-dV_DSM_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));
[output_1] = NL_interpolator_of( r_GA1 , rA1 , v_dep_GA1 , vA1 , N_rev1 , TOF1 , Mass_at_beginning_LT ,sim.PS.Isp , sim );
%     if ~isnan(output_1.Thrust(1,1)) && abs(output_1.t(end) - TOF1) < tol_TOF % if is not nan -> it's a valid solution
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    %penalty_T_leg1 = abs(max(output_1.T_magn)) - sim.max_Available_Thrust;
    penalty_T_leg1 = max(abs(output_1.T_magn)) - sim.max_Available_Thrust;
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
v_arr_GA1_dim = v_arr_GA1.*sim.DU./sim.TU;
v_dep_GA1_dim = v_dep_GA1.*sim.DU./sim.TU;
[delta_v_p_1,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGA1_dim, v_arr_GA1_dim, v_dep_GA1_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_1), 'Not found')
    penalty_dv_GA1 = 10;
end

% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end) - sim.M_pods; % yes pods here, released on Ast1
[output_2] = NL_interpolator_of( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
%     if ~isnan(output_2.Thrust(1,1)) && abs(output_2.t(end) - TOF2) < tol_TOF  % if is not nan -> it's a valid solution
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = max(abs(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

%% SC 2
% ballistic stuff before GA
% DV calculation with lambert
% needs km, km/s, km^3/s^2; it has also h as 7th element
% output km, rad, km^2/s
coe_transf_orbit_7_el = coe_from_sv(r_EA*sim.DU,V_EA_DSM_dep*sim.DU/sim.TU,ksun); 
coe_transf_orbit = coe_transf_orbit_7_el(1:6);
% determination of point where happens the 2nd DSM
DSM2_theta_in_transfer_orbit = coe_transf_orbit(6)+thetaD2;
coe_DSM2 = [coe_transf_orbit(1:5),DSM2_theta_in_transfer_orbit];
[rDSM2, vDSM2] = sv_from_coe(coe_DSM2,ksun); % output; km, km/s
rDSM2 = rDSM2/sim.DU;
vDSM2 = vDSM2/sim.DU*sim.TU;

% -------------------- times of stuff happening for SC2 ----------------- %
% ------ possible only after the common transfer orbit has been defined ------ %
% -- DSM2
% 2nd spacecraft characteristic times
% TOF_DSM2 = x(16); % It's the time from Earth to DSM
% ksun    Gravity constant of the Sun [km^3/s^2]
TOF_DSM2_in_seconds = time_law(coe_transf_orbit(6),DSM2_theta_in_transfer_orbit,coe_transf_orbit(1),coe_transf_orbit(2),ksun);
TOF_DSM2 = TOF_DSM2_in_seconds/sim.TU;
TOFGA_a = x(17); % DSM to GA
MJDPGA_a = MJD01 + TOF_DSM2 + TOFGA_a; % passage time GA
TOFa = x(18); % tof sc1 GA to Ast1 
MJDAa = MJDPGA_a + TOFa; % mjd2000 arrive of 1st sc on ast 1
CTa = x(25);
MJDDa = MJDAa + CTa; % departure from 1st asteroid
TOFb = x(26); % tof sc1 from ast 1 to 2nd asteroid
MJDAb = MJDDa + TOFb; % mjd2000 passage of SC1 on ast 2

% GA SC2
MJDPGAa_dim = MJDPGA_a*sim.TU/(3600*24);
[kep_GAa,ksun] = uplanet(MJDPGAa_dim, sim.ID_FLYBY);
[r_GAa, v_GAa] = sv_from_coe(kep_GAa,ksun);
r_GAa = r_GAa/sim.DU;
v_GAa = v_GAa/sim.DU*sim.TU;

% ARRIVAL at a_th ast
MJDAa_dim = MJDAa*sim.TU/(3600*24);
[kep_ast_Aa] = uNEO2(MJDAa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rAa, vAa] = sv_from_coe(kep_ast_Aa,ksun); % km, km/s
rAa = rAa/sim.DU;
vAa = vAa/sim.DU*sim.TU;
% DEPARTURE at a_th ast
MJDDa_dim = MJDDa*sim.TU/(3600*24);
[kep_ast_Da] = uNEO2(MJDDa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rDa, vDa] = sv_from_coe(kep_ast_Da,ksun); % km, km/s
rDa = rDa/sim.DU;
vDa = vDa/sim.DU*sim.TU;

% passage at b_th ast
MJDAb_dim = MJDAb*sim.TU/(3600*24);
[kep_ast_Ab] = uNEO2(MJDAb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;

% Gravity assist SC2
v_arr_GAa = v_GAa + v_inf_magn_GAa*[cos(el_in_GAa)*cos(az_in_GAa); cos(el_in_GAa)*sin(az_in_GAa); sin(el_in_GAa)];

% -------- manoeuvres after DSM for SC2 --------------- %

% DSM -> GA
[~,~,~,~,V_DSM2_GA_dep,V_DSM2_GA_arr,~,~] = lambertMR(rDSM2,r_GAa,TOFGA_a,sim.mu,0,0,0,0);
V_DSM2_GA_dep = V_DSM2_GA_dep';
V_DSM2_GA_arr = V_DSM2_GA_arr';
dV_DSM2_GA = sqrt((V_DSM2_GA_dep(1)-vDSM2(1))^2+(V_DSM2_GA_dep(2)-vDSM2(2))^2+(V_DSM2_GA_dep(3)-vDSM2(3))^2);
if dV_DSM2_GA < sim.dV_DSM_max % vinf that the launcher can give max 
    pen_dV_DSM2 = 0;
else
    pen_dV_DSM2 = 100*abs(dV_DSM2_GA - sim.dV_DSM_max); % penalty like, but not discard a priori
end

% 1st leg - GA -> Ast 1
v_dep_GAa = v_GAa + v_inf_magn_GAa*[cos(el_out_GAa)*cos(az_out_GAa); cos(el_out_GAa)*sin(az_out_GAa); sin(el_out_GAa)];
Mass_at_beginning_LT2 = sim.M2*exp(-dV_DSM2_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));
[output_a] = NL_interpolator_of( r_GAa , rAa , v_dep_GAa , vAa , N_reva , TOFa , Mass_at_beginning_LT2 ,sim.PS.Isp , sim );
%     if ~isnan(output_1.Thrust(1,1)) && abs(output_1.t(end) - TOF1) < tol_TOF % if is not nan -> it's a valid solution
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
    %penalty_T_leg1 = abs(max(output_1.T_magn)) - sim.max_Available_Thrust;
    penalty_T_lega = max(abs(output_a.T_magn)) - sim.max_Available_Thrust;
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
v_arr_GAa_dim = v_arr_GAa.*sim.DU./sim.TU;
v_dep_GAa_dim = v_dep_GAa.*sim.DU./sim.TU;
[delta_v_p_a,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGAa_dim, v_arr_GAa_dim, v_dep_GAa_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_a), 'Not found')
    penalty_dv_GAa = 10;
end

% bth leg - Asta -> Astb
M_start_bth_leg = output_a.m(end) - sim.M_pods; % yes pods here, released on Ast1
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , M_start_bth_leg ,sim.PS.Isp , sim );
%     if ~isnan(output_2.Thrust(1,1)) && abs(output_2.t(end) - TOF2) < tol_TOF  % if is not nan -> it's a valid solution
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = max(abs(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

%% masses and obj function
Isp = sim.PS.Isp*sim.TU;

% SC1
mass_depleted_DSM1 = sim.M1 - Mass_at_beginning_LT;
mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);

tot_mass_depleted_SC1 = mass_depleted_DSM1+mass_depleted_Leg1+mass_depleted_Leg2;
mass_dry_and_pods_SC1 = sim.M1 - tot_mass_depleted_SC1;
mass_fract_SC1 = (sim.M1 - mass_dry_and_pods_SC1)/sim.M1;

sol.dV_DSM1 = dV_DSM_GA*sim.DU/sim.TU*1e3;
sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

% gravity assist quantities
sol.delta_V_p1 = delta_v_p_1;
sol.dV_gained_GA1 = (sqrt((v_dep_GA1(1)-V_DSM_GA_arr(1))^2+(v_dep_GA1(2)-V_DSM_GA_arr(2))^2+(v_dep_GA1(3)-V_DSM_GA_arr(3))^2) - sol.delta_V_p1)*sim.DU/sim.TU;

% SC2
mass_depleted_DSM2 = sim.M2 - Mass_at_beginning_LT2;
mass_depleted_Lega = output_a.m(1) - output_a.m(end);
mass_depleted_Legb = output_b.m(1) - output_b.m(end);

tot_mass_depleted_SC2 = mass_depleted_DSM2+mass_depleted_Lega+mass_depleted_Legb;
mass_dry_and_pods_SC2 = sim.M2 - tot_mass_depleted_SC2;
mass_fract_SC2 = (sim.M2 - mass_dry_and_pods_SC2)/sim.M2;

sol.dV_DSM2 = dV_DSM2_GA*sim.DU/sim.TU*1e3;
sol.dV_associated_Lega = -g0*Isp*log(output_a.m(end)/output_a.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_Legb = -g0*Isp*log(output_b.m(end)/output_b.m(1)); % -ve*ln(m_final/m_initial)

% gravity assist quantities
sol.delta_V_pa = delta_v_p_a;
sol.dV_gained_GAa = (sqrt((v_dep_GAa(1)-V_DSM2_GA_arr(1))^2+(v_dep_GAa(2)-V_DSM2_GA_arr(2))^2+(v_dep_GAa(3)-V_DSM2_GA_arr(3))^2) - sol.delta_V_pa)*sim.DU/sim.TU;

penalty_MF_unfeasible = 0;
if mass_fract_SC1 <= 0 || mass_fract_SC1 >= 1 || mass_fract_SC2 <= 0 || mass_fract_SC2 >= 1
    penalty_MF_unfeasible = 100;
end

avg_mass_fraction = (mass_fract_SC1+mass_fract_SC2)/2;
MF = max(mass_fract_SC1,mass_fract_SC2) + abs(mass_fract_SC1 - avg_mass_fraction) + ...
        abs(mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali

sol.obj_fun = MF + penalty_MAX_DURATION + penalty_MF_unfeasible + ...
    20*(dV_DSM_GA + dV_DSM2_GA) + penalty_dV_C3 + ...
    penalty_dv_GA1 + penalty_dv_GAa + ...
    20*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb + ....
    pen_dV_DSM + pen_dV_DSM2;
    
%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.DSM1 = rDSM1;
r_encounter.GA1 = r_GA1;
r_encounter.astA1 = rA1;
r_encounter.astD1 = rD1;
r_encounter.astA2 = rA2;
r_encounter.DSM2 = rDSM2;
r_encounter.GAa = r_GAa;
r_encounter.astAa = rAa;
r_encounter.astDa = rDa;
r_encounter.astAb = rAb;

v_encounter.EA = v_EA;
% v_encounter.DSM1 = vDSM1;
v_encounter.GA1 = v_GA1;
v_encounter.astA1 = vA1;
v_encounter.astD1 = vD1;
v_encounter.astA2 = vA2;
% v_encounter.DSM2 = vDSM2;
v_encounter.GAa = v_GAa;
v_encounter.astAa = vAa;
v_encounter.astDa = vDa;
v_encounter.astAb = vAb;

%% Porcherie
% t_span_EA_DSM1 = linspace(0,TOF_DSM1*sim.TU,sim.n_sol);
% t_span_DSM1_GA1 = linspace(TOF_DSM1*sim.TU,(TOF_DSM1+TOFGA_1)*sim.TU,sim.n_sol);
% t_span_CT1 = linspace((TOF_DSM1+TOFGA_1)*sim.TU+output_1.t(end),(TOF_DSM1+TOFGA_1)*sim.TU+output_1.t(end)+CT1,sim.n_sol);
% t_span_EA_DSM2 = linspace(0,TOF_DSM2*sim.TU,sim.n_sol);
% t_span_DSM2_GA2 = linspace(TOF_DSM2*sim.TU,(TOF_DSM2+TOFGA_a)*sim.TU,sim.n_sol);
% t_span_CTa = linspace((TOF_DSM2+TOFGA_a)*sim.TU+output_a.t(end),(TOF_DSM2+TOFGA_a)*sim.TU+output_a.t(end)+CTa,sim.n_sol);
% only LT plotting
t_span_CT1 = linspace(output_1.t(end),output_1.t(end)+CT1,sim.n_sol);
t_span_CTa = linspace(output_a.t(end),output_a.t(end)+CTa,sim.n_sol);
mCT1 = ones(sim.n_sol,1).*output_1.m(end);
mCTa = ones(sim.n_sol,1).*output_a.m(end);
TCT1 = zeros(sim.n_sol,3); 
TCTa = zeros(sim.n_sol,3);

output.t_SC1            = [output_1.t; t_span_CT1'; t_span_CT1(end)+output_2.t];
output.t_SC2            = [output_a.t; t_span_CTa'; t_span_CTa(end)+output_b.t];
output.m_SC1            = [output_1.m; mCT1; output_2.m];
output.m_SC2            = [output_a.m; mCTa; output_b.m];
output.Thrust_SC1       = [output_1.Thrust; TCT1; output_2.Thrust];
output.T_magn_SC1       = sqrt(output.Thrust_SC1(:,1).^2 + output.Thrust_SC1(:,3).^2);
output.Thrust_SC2       = [output_a.Thrust; TCTa;output_b.Thrust];
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