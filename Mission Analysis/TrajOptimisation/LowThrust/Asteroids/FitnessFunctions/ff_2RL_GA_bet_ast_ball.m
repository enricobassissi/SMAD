function obj_fun = ff_2RL_GA_bet_ast_ball(x,sim,data,power_propulsion_data)
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
     (16) CT2,
% GA Stuff SC1
TOFGA (17)
azimuth out (POST-GRAVITY ASSIST)(22)
elevation (POST-GRAVITY ASSIST)(23)
% GA Stuff SC2
TOFGA (24)
azimuth out (POST-GRAVITY ASSIST)(29)
elevation (POST-GRAVITY ASSIST)(30)]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDA1 = MJD01 + TOF1; % mjd2000 arrive of 1st sc on ast 1
CT1 = x(15);
MJDD1 = MJDA1 + CT1; % departure from 1st asteroid
TOFGA_1 = x(17);
MJDPGA_1 = MJDD1 + TOFGA_1; % passage time GA
TOF2 = x(3); % tof sc1 to 2nd asteroid
MJDA2 = MJDPGA_1 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = x(4); % tof sc2 to 1st asteroid
MJDAa = MJD01 + TOFa; 
CTa = x(16);
MJDDa = MJDAa + CTa;
TOFGA_a = x(20);
MJDPGA_a = MJDDa + TOFGA_a; % passage time GA
TOFb =  x(5); % tof sc2 to 2nd asteroid
MJDAb = MJDPGA_a + TOFb; 

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

% GA stuff SC1
% the condition IN is determined by the Lambert Leg incoming
% we change the late
% OUT
az_out_GA1 = x(18);
el_out_GA1 = x(19);

% GA stuff SC2
% OUT
az_out_GAa = x(21);
el_out_GAa = x(22);

%% choosing which asteroid to visit
% 1ST SPACECRAFT ASTEROID OBJECTIVES
IDP1 = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% % 2ND SPACECRAFT ASTEROID OBJECTIVES
IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
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

% GA SC1
MJDPGA1_dim = MJDPGA_1*sim.TU/(3600*24);
[kep_GA1,ksun] = uplanet(MJDPGA1_dim, sim.ID_FLYBY);
[r_GA1, v_GA1] = sv_from_coe(kep_GA1,ksun);
r_GA1 = r_GA1/sim.DU;
% v_GA1 = v_GA1/sim.DU*sim.TU;

% ARRIVAL at 2nd ast
MJDA2_dim = MJDA2*sim.TU/(3600*24);
[kep_ast_A2] = uNEO2(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;

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

% GA SC2
MJDPGAa_dim = MJDPGA_a*sim.TU/(3600*24);
[kep_GAa,ksun] = uplanet(MJDPGAa_dim, sim.ID_FLYBY);
[r_GAa, v_GAa] = sv_from_coe(kep_GAa,ksun);
r_GAa = r_GAa/sim.DU;
% v_GAa = v_GAa/sim.DU*sim.TU;

% passage at b_th ast
MJDAb_dim = MJDAb*sim.TU/(3600*24);
[kep_ast_Ab] = uNEO2(MJDAb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;

%% Launcher departure variable
v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

%% NLI %%
tol_TOF = 1; % 1 TU means approx 60 days
penalty_dv_GA1 = 0; penalty_dv_GAa = 0;
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;

Isp_secondary_prop = 255; % ADN
g0 = 9.81;

%% SC1 
% 1st leg - Earth -> Ast 1
[output_1] = NL_interpolator_of( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , sim.M1 ,sim.PS.Isp , sim );
% if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
%     penalty_T_leg1 = max(abs(output_1.T_magn)) - sim.max_Available_Thrust;
% end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end
if ~isnan(output_1.Href)
    % t in years
    t = output_1.t*sim.TU/(3600*24*365);
    % trajectory heliocentric
    r3  = [output_1.r.*cos(output_1.theta) output_1.r.*sin(output_1.theta) output_1.z];
    R3 = rotate_local2ecplitic(r_EA,r3,sim.n_sol,output_1.Href);
    % thrust in heliocentric
    Tlocal_transf_orbit  = [-output_1.Thrust(:,1).*sin(output_1.theta), ...
        output_1.Thrust(:,1).*cos(output_1.theta), output_1.Thrust(:,3)];
    Thrust_Helio = rotate_local2ecplitic(r_EA,Tlocal_transf_orbit,sim.n_sol,output_1.Href);
    T_magn_Helio = vecnorm(Thrust_Helio,2,2).*1e3; % mN
    %thrust available
    T_available = available_thrust(t, R3, Thrust_Helio, power_propulsion_data); % mN
    if max(T_magn_Helio - T_available) > 0
        penalty_T_leg1 = max(T_magn_Helio - T_available);
    end
else
    penalty_T_leg1 = 100;
end

% 2nd leg - Ast1 -> GA SC1
[~,~,~,~,V_D1_GA_dep,V_D1_GA_arr,~,~] = lambertMR(rD1,r_GA1,TOFGA_1,sim.mu,0,0,0,0);
V_D1_GA_dep = V_D1_GA_dep';
V_D1_GA_arr = V_D1_GA_arr';
dV_D1_GA = sqrt((V_D1_GA_dep(1)-vD1(1))^2+(V_D1_GA_dep(2)-vD1(2))^2+(V_D1_GA_dep(3)-vD1(3))^2);
if dV_D1_GA < sim.dV_man_max % vinf that the launcher can give max 
    pen_dV_man1 = 0;
else
    pen_dV_man1 = 1e4*abs(dV_D1_GA - sim.dV_man_max); % penalty like, but not discard a priori
end

Mass_before_burn_to_GA1 = (output_1.m(end) - sim.M_pods); % yes pods here, released on Ast1
Mass_GA1_leg = Mass_before_burn_to_GA1*exp(-dV_D1_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop)); % depletion of mass to start going to GA
 % [output_GA1] = NL_interpolator_of( rD1 , r_GA1 , vD1 , v_arr_GA1 , N_rev_GA1 , TOFGA_1 , M_start_GA1_leg ,sim.PS.Isp , sim );
% if max(abs(output_GA1.T_magn)) > sim.max_Available_Thrust
%     %penalty_T_leg2 = abs(max(output_2.T_magn)) - sim.max_Available_Thrust;
%     penalty_T_legGA1 = max(abs(output_GA1.T_magn)) - sim.max_Available_Thrust;
% end
% if abs(output_GA1.t(end) - TOFGA_1) > tol_TOF
%     penalty_TOF_legGA1 = abs(output_GA1.t(end) - TOFGA_1);
% end

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
v_arr_GA1_dim = V_D1_GA_arr.*sim.DU./sim.TU;
v_dep_GA1 = norm(V_D1_GA_arr)*[cos(el_out_GA1)*cos(az_out_GA1); cos(el_out_GA1)*sin(az_out_GA1); sin(el_out_GA1)];
v_dep_GA1_dim = v_dep_GA1.*sim.DU./sim.TU;
[delta_v_p_1,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGA1_dim, v_arr_GA1_dim, v_dep_GA1_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_1), 'Not found')
    penalty_dv_GA1 = 100;
end

% 3rd leg - GA SC1 -> Ast 2
[output_2] = NL_interpolator_of( r_GA1 , rA2 , v_dep_GA1 , vA2 , N_rev2 , TOF2 , Mass_GA1_leg ,sim.PS.Isp , sim );
% if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
%     penalty_T_leg2 = max(abs(output_2.T_magn)) - sim.max_Available_Thrust;
% end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end
if ~isnan(output_2.Href)
    % t in years
    t = output_2.t*sim.TU/(3600*24*365);
    % trajectory heliocentric
    r3  = [output_2.r.*cos(output_2.theta) output_2.r.*sin(output_2.theta) output_2.z];
    R3 = rotate_local2ecplitic(r_GA1,r3,sim.n_sol,output_2.Href);
    % thrust in heliocentric
    Tlocal_transf_orbit  = [-output_2.Thrust(:,1).*sin(output_2.theta), ...
        output_2.Thrust(:,1).*cos(output_2.theta), output_2.Thrust(:,3)];
    Thrust_Helio = rotate_local2ecplitic(r_GA1,Tlocal_transf_orbit,sim.n_sol,output_2.Href);
    T_magn_Helio = vecnorm(Thrust_Helio,2,2).*1e3; % mN
    %thrust available
    T_available = available_thrust(t, R3, Thrust_Helio, power_propulsion_data); % mN
    if max(T_magn_Helio - T_available) > 0
        penalty_T_leg2 = max(T_magn_Helio - T_available);
    end
else
    penalty_T_leg2 = 100;
end

%% SC 2
% a_th leg - Earth -> Ast_a
[output_a] = NL_interpolator_of( r_EA , rAa , v_dep , vAa , N_reva , TOFa , sim.M2 ,sim.PS.Isp , sim );
% if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
%     penalty_T_lega = max(abs(output_a.T_magn)) - sim.max_Available_Thrust;
% end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end
if ~isnan(output_a.Href)
    % t in years
    t = output_a.t*sim.TU/(3600*24*365);
    % trajectory heliocentric
    r3  = [output_a.r.*cos(output_a.theta) output_a.r.*sin(output_a.theta) output_a.z];
    R3 = rotate_local2ecplitic(r_EA,r3,sim.n_sol,output_a.Href);
    % thrust in heliocentric
    Tlocal_transf_orbit  = [-output_a.Thrust(:,1).*sin(output_a.theta), ...
        output_a.Thrust(:,1).*cos(output_a.theta), output_a.Thrust(:,3)];
    Thrust_Helio = rotate_local2ecplitic(r_EA,Tlocal_transf_orbit,sim.n_sol,output_a.Href);
    T_magn_Helio = vecnorm(Thrust_Helio,2,2).*1e3; % mN
    %thrust available
    T_available = available_thrust(t, R3, Thrust_Helio, power_propulsion_data); % mN
    if max(T_magn_Helio - T_available) > 0
        penalty_T_lega = max(T_magn_Helio - T_available);
    end
else
    penalty_T_lega = 100;
end

% GAa Leg - Ast_a -> GAa SC 2
[~,~,~,~,V_Da_GA_dep,V_Da_GA_arr,~,~] = lambertMR(rDa,r_GAa,TOFGA_a,sim.mu,0,0,0,0);
V_Da_GA_dep = V_Da_GA_dep';
V_Da_GA_arr = V_Da_GA_arr';
dV_Da_GA = sqrt((V_Da_GA_dep(1)-vDa(1))^2+(V_Da_GA_dep(2)-vDa(2))^2+(V_Da_GA_dep(3)-vDa(3))^2);
if dV_Da_GA < sim.dV_man_max % vinf that the launcher can give max 
    pen_dV_mana = 0;
else
    pen_dV_mana = 1e4*abs(dV_Da_GA - sim.dV_man_max); % penalty like, but not discard a priori
end

Mass_before_burn_to_GAa = (output_a.m(end) - sim.M_pods); % yes pods here, released on Ast1
Mass_GAa_leg = Mass_before_burn_to_GAa*exp(-dV_Da_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop)); % depletion of mass to start going to GA

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
v_arr_GAa_dim = V_Da_GA_arr.*sim.DU./sim.TU;
v_dep_GAa = norm(V_Da_GA_arr)*[cos(el_out_GAa)*cos(az_out_GAa); cos(el_out_GAa)*sin(az_out_GAa); sin(el_out_GAa)];
v_dep_GAa_dim = v_dep_GAa.*sim.DU./sim.TU;
[delta_v_p_a,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGAa_dim, v_arr_GAa_dim, v_dep_GAa_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_a), 'Not found')
    penalty_dv_GAa = 100;
end

% b_th Leg - GAa -> Ast_b
[output_b] = NL_interpolator_of( r_GAa , rAb , vDa , vAb , N_revb , TOFb , Mass_GAa_leg ,sim.PS.Isp , sim );
% if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
%     penalty_T_legb = max(abs(output_b.T_magn)) - sim.max_Available_Thrust;
% end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end
if ~isnan(output_b.Href)
    % t in years
    t = output_b.t*sim.TU/(3600*24*365);
    % trajectory heliocentric
    r3  = [output_b.r.*cos(output_b.theta) output_b.r.*sin(output_b.theta) output_b.z];
    R3 = rotate_local2ecplitic(r_GAa,r3,sim.n_sol,output_b.Href);
    % thrust in heliocentric
    Tlocal_transf_orbit  = [-output_b.Thrust(:,1).*sin(output_b.theta), ...
        output_b.Thrust(:,1).*cos(output_b.theta), output_b.Thrust(:,3)];
    Thrust_Helio = rotate_local2ecplitic(r_GAa,Tlocal_transf_orbit,sim.n_sol,output_b.Href);
    T_magn_Helio = vecnorm(Thrust_Helio,2,2).*1e3; % mN
    %thrust available
    T_available = available_thrust(t, R3, Thrust_Helio, power_propulsion_data); % mN
    if max(T_magn_Helio - T_available) > 0
        penalty_T_legb = max(T_magn_Helio - T_available);
    end
else
    penalty_T_legb = 100;
end

%% masses and obj function
% SC1
mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
% mass_depleted_GA1_inj = Mass_before_burn_to_GA1 - Mass_GA1_leg;
mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);

% tot_mass_depleted_SC1 = mass_depleted_Leg1+mass_depleted_GA1_inj+mass_depleted_Leg2;
LT_mass_depleted_SC1 = mass_depleted_Leg1+mass_depleted_Leg2;
% mass_dry_and_pods_SC1 = output_1.m(1) - tot_mass_depleted_SC1;
% mass_fract_SC1 = (output_1.m(1) - mass_dry_and_pods_SC1)/output_1.m(1);
mass_fract_SC1 = LT_mass_depleted_SC1/output_1.m(1);

% SC2
mass_depleted_Lega = output_a.m(1) - output_a.m(end);
% mass_depleted_GAa_inj = Mass_before_burn_to_GAa - Mass_GAa_leg;
mass_depleted_Legb = output_b.m(1) - output_b.m(end);

% tot_mass_depleted_SC2 = mass_depleted_Lega+mass_depleted_GAa_inj+mass_depleted_Legb;
% mass_dry_and_pods_SC2 = output_a.m(1) - tot_mass_depleted_SC2;
LT_mass_depleted_SC2 = mass_depleted_Lega+mass_depleted_Legb;
% mass_fract_SC2 = (output_a.m(1) - mass_dry_and_pods_SC2)/output_a.m(1);
mass_fract_SC2 = LT_mass_depleted_SC2/output_a.m(1);

%     if mass_fract_SC1 > 0 && mass_fract_SC1 < 1 
%         if mass_fract_SC2 > 0 && mass_fract_SC2 < 1 
avg_mass_fraction = (mass_fract_SC1+mass_fract_SC2)/2;
MF = max(mass_fract_SC1,mass_fract_SC2) + abs(mass_fract_SC1 - avg_mass_fraction) + ...
        abs(mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali

obj_fun = MF + penalty_MAX_DURATION + penalty_dv_GA1 + penalty_dv_GAa + ...
    20*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb) + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb + ...
    pen_dV_man1 + pen_dV_mana;

end