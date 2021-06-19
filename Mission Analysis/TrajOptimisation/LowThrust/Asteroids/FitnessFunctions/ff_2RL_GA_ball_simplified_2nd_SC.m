function obj_fun = ff_2RL_GA_ball_simplified_2nd_SC(x,sim,data,power_propulsion_data)
%{ 

input variable vector
% x = [
% thetaD2 (1)
% TOFGAa (2), from DSM2 to GA
% azimuth out (POST-GRAVITY ASSIST)(3)
% elevation (POST-GRAVITY ASSIST)(4)
% TOFa (5)
% N REV a (6)
% Coasting time a (7)
% TOFb (8)
% N REV b (9)

%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% -- DSM SC2
thetaD2_rel = x(1);

% setting the input times
MJD01 = sim.departure_adim; % departure time for both the sc
TOFGA_a = x(2);
TOFa = x(5);
CTa = x(7);
TOFb = x(8);

% N REVa
N_reva = x(6);
% N REVb
N_revb = x(9);

% GA stuff SC2
% OUT
az_out_GAa = x(3);
el_out_GAa = x(4);

%% choosing which asteroid to visit
IDP2 = x(10);
asteroid_a = data.PermutationMatrix_2SC(IDP2,1);
asteroid_b = data.PermutationMatrix_2SC(IDP2,2);

% arbitrarly fixed, they're nice in sequence
% asteroid_a = "2010UJ";
% asteroid_b = "2016CK137";

%% FUEL CONSUMPTIONS %%
tol_TOF = 1; % 1 TU means approx 60 days
penalty_dv_GAa = 0;
penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_lega = 0; penalty_TOF_legb = 0;

Isp_secondary_prop = 255; % ADN
g0 = 9.81;

%% SC 2
% ballistic stuff before GA
% DV calculation with lambert
% needs km, km/s, km^3/s^2; it has also h as 7th element
% output km, rad, km^2/s
coe_transf_orbit_7_el = coe_from_sv(sim.r_EA_dep*sim.DU,sim.V_EA_DSM_dep*sim.DU/sim.TU,sim.mu_dim); 
coe_transf_orbit = coe_transf_orbit_7_el(1:6); % coe of the orbit of transfer from EA departure to DSM1, but at r_EA moment!
% determination of point where happens the 2nd DSM -> at moment of DSM1 +/- something (about 10°)
DSM2_theta_in_transfer_orbit = coe_transf_orbit(6)+sim.theta_DSM1+thetaD2_rel;
coe_DSM2 = [coe_transf_orbit(1:5),DSM2_theta_in_transfer_orbit];
[rDSM2, vDSM2] = sv_from_coe(coe_DSM2,sim.mu_dim); % input: km, rad; output; km, km/s
rDSM2 = rDSM2/sim.DU;
vDSM2 = vDSM2/sim.DU*sim.TU;

% -------------------- times of stuff happening for SC2 ----------------- %
% ------ possible only after the common transfer orbit has been defined ------ %
% -- DSM2
% 2nd spacecraft characteristic times
% ksun    Gravity constant of the Sun [km^3/s^2]
penalty_dsm_hyp = 0;
if coe_transf_orbit(1) > 0
    TOF_DSM2_in_seconds = time_law(coe_transf_orbit(6),DSM2_theta_in_transfer_orbit,coe_transf_orbit(1),coe_transf_orbit(2),sim.mu_dim);
    TOF_DSM2 = TOF_DSM2_in_seconds/sim.TU;
else
    TOF_DSM2 = 10;
    penalty_dsm_hyp = 100;
end

MJDPGA_a = MJD01 + TOF_DSM2 + TOFGA_a; % passage time GA
MJDAa = MJDPGA_a + TOFa; % mjd2000 arrive of 1st sc on ast 1
MJDDa = MJDAa + CTa; % departure from 1st asteroid
MJDAb = MJDDa + TOFb; % mjd2000 passage of SC1 on ast 2

% GA SC2
MJDPGAa_dim = MJDPGA_a*sim.TU/(3600*24);
[kep_GAa,ksun] = uplanet(MJDPGAa_dim, 3);
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

% -------- manoeuvres after DSM for SC2 --------------- %
% DSM -> GA
[~,~,~,~,V_DSM2_GA_dep,V_DSM2_GA_arr,~,~] = lambertMR(rDSM2,r_GAa,TOFGA_a,sim.mu,0,0,0,0);
V_DSM2_GA_dep = V_DSM2_GA_dep';
V_DSM2_GA_arr = V_DSM2_GA_arr';
dV_DSM2_GA = sqrt((V_DSM2_GA_dep(1)-vDSM2(1))^2+(V_DSM2_GA_dep(2)-vDSM2(2))^2+(V_DSM2_GA_dep(3)-vDSM2(3))^2);
if dV_DSM2_GA < sim.dV_DSM_max % vinf that the launcher can give max 
    pen_dV_DSM2 = 0;
else
    pen_dV_DSM2 = 1000*abs(dV_DSM2_GA - sim.dV_DSM_max); % penalty like, but not discard a priori
end

% ath leg - GA -> Ast a
v_inf_minus = V_DSM2_GA_arr - v_GAa;
v_inf_plus = norm(v_inf_minus)*[cos(el_out_GAa)*cos(az_out_GAa); cos(el_out_GAa)*sin(az_out_GAa); sin(el_out_GAa)];
v_dep_GAa = v_inf_plus + v_GAa; 
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
v_arr_GAa_dim = V_DSM2_GA_arr.*sim.DU./sim.TU;
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
% SC2
mass_depleted_DSM2 = sim.M2 - Mass_at_beginning_LT2;
mass_depleted_Lega = output_a.m(1) - output_a.m(end);
mass_depleted_Legb = output_b.m(1) - output_b.m(end);

tot_mass_depleted_SC2 = mass_depleted_DSM2+mass_depleted_Lega+mass_depleted_Legb;
mass_dry_and_pods_SC2 = sim.M2 - tot_mass_depleted_SC2;
mass_fract_SC2 = (sim.M2 - mass_dry_and_pods_SC2)/sim.M2;

penalty_MF_unfeasible = 0;
if mass_fract_SC2 <= 0 || mass_fract_SC2 >= 1
    penalty_MF_unfeasible = 100;
end

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if TOF_DSM2+TOFGA_a+TOFa+CTa+TOFb > max_duration
    penalty_MAX_DURATION = TOF_DSM2+TOFGA_a+TOFa+CTa+TOFb - max_duration; % 12 years max mission time 
end

obj_fun = mass_fract_SC2 + penalty_MAX_DURATION + penalty_MF_unfeasible + ...
    20*dV_DSM2_GA + ...
    penalty_dv_GAa + ...
    20*(penalty_T_lega + penalty_T_legb) + ...
    penalty_TOF_lega + penalty_TOF_legb + ...
    pen_dV_DSM2 + penalty_dsm_hyp;

end

