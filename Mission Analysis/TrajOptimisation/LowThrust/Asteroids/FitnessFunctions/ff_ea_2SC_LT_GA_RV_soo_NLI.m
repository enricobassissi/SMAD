function obj_fun = ff_ea_2SC_LT_GA_RV_soo_NLI(x,sim,data)
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
    TOFGA (17) %GA SC1
    N REV GA (18)
    v_inf_magn GA (19)
    azimuth GAin (20)
    elevation GAin (21)
    azimuth GAout (22)
    elevation GAout (23),
    TOFGA2 (24) % GA SC2
    N REV GA2 (25)
    v_inf_magn GA2 (26)
    azimuth GA2in (27)
    elevation GA2in (28)
    azimuth GA2out (29)
    elevation GA2out (30)]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOFGA_SC1 = x(17);
MJDGA1 = MJD01 + TOFGA_SC1; % mjd2000 of gravity assist
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDA1 = MJDGA1 + TOF1; % mjd2000 arrival of 1st sc on ast 1
CT1 = x(15);
MJDD1 = MJDA1 + CT1;
TOF2 = x(3); % tof sc1 to 2nd asteroid
MJDA2 = MJDD1 + TOF2; % mjd2000 arrival of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFGA_SC2 = x(24);
MJDGA2 = MJD01 + TOFGA_SC2; % mjd2000 of gravity assist
TOFa = x(4); % tof sc2 to 1st asteroid
MJDAa = MJDGA2 + TOFa; 
CTa = x(16);
MJDDa = MJDAa + CTa;
TOFb =  x(5); % tof sc2 to 2nd asteroid
MJDAb = MJDDa + TOFb; 

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if max(TOFGA_SC1+TOF1+CT1+TOF2,TOFGA_SC2+TOFa+CTa+TOFb) > max_duration
    penalty_MAX_DURATION = max(TOFGA_SC1+TOF1+CT1+TOF2,TOFGA_SC2+TOFa+CTa+TOFb) - max_duration; % 12 years max mission time 
end

% N REV1
N_rev1 = x(6);
% N REV2
N_rev2 = x(7);
% N REVa
N_reva = x(8);
% N REVb
N_revb = x(9);
% N REVGA1
N_revGA1 = x(18);
% N REVGA2
N_revGA2 = x(25);

% C3 launcher
v_inf_magn = x(12);
az = x(13);
elev = x(14);

%% Gravity assists stuff
% SC1
% GA stuff IN
v_inf_magnGA = x(19);
az_GAin = x(20);
el_GAin = x(21);
% GA stuff out
az_GAout = x(22);
el_GAout = x(23);
% SC2
% GA stuff IN
v_inf_magnGA2 = x(26);
az_GA2in = x(27);
el_GA2in = x(28);
% GA stuff out
az_GA2out = x(29);
el_GA2out = x(30);

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

%% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% GA SC1
MJDGA1_dim = MJDGA1*sim.TU/(3600*24);
[kep_GA1,ksun] = uplanet(MJDGA1_dim, sim.ID_FLYBY);
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

% GA SC2
MJDGA2_dim = MJDGA2*sim.TU/(3600*24);
[kep_GA2,ksun] = uplanet(MJDGA2_dim, sim.ID_FLYBY);
[r_GA2, v_GA2] = sv_from_coe(kep_GA2,ksun);
r_GA2 = r_GA2/sim.DU;
v_GA2 = v_GA2/sim.DU*sim.TU;

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

%% Launcher departure variable
%     z = r .* sin(elev);
%     rcoselev = r .* cos(elev);
%     x = rcoselev .* cos(az);
%     y = rcoselev .* sin(az);
v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

% Gravity assist
v_arr_GA1 = v_GA1 + v_inf_magnGA*[cos(el_GAin)*cos(az_GAin); cos(el_GAin)*sin(az_GAin); sin(el_GAin)];
v_arr_GA2 = v_GA2 + v_inf_magnGA2*[cos(el_GA2in)*cos(az_GA2in); cos(el_GA2in)*sin(az_GA2in); sin(el_GA2in)];

%% NLI
% intro
tol_TOF = 1; % 1 TU means approx 60 days
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;
penalty_T_legGA1 = 0; penalty_TOF_legGA1 = 0; penalty_dv_GA1 = 0; 
penalty_T_legGA2 = 0; penalty_TOF_legGA2 = 0; penalty_dv_GA2 = 0; 

% SC1 
% 0th leg - Earth -> GA1
[output_GA1] = NL_interpolator_of( r_EA , r_GA1 , v_dep , v_arr_GA1 , N_revGA1 , TOFGA_SC1 , sim.M1 ,sim.PS.Isp , sim );
if abs(max(output_GA1.T_magn)) > sim.max_Available_Thrust
    penalty_T_legGA1 = 10*abs(max(output_GA1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_GA1.t(end) - TOFGA_SC1) > tol_TOF
    penalty_TOF_legGA1 = abs(output_GA1.t(end) - TOFGA_SC1);
end
% gravity assist
v_dep_GA1 = v_GA1 + v_inf_magnGA*[cos(el_GAout)*cos(az_GAout); cos(el_GAout)*sin(az_GAout); sin(el_GAout)];
if sim.ID_FLYBY == 3
    RPlanet_flyby = astroConstants(23); % Radius_Earth, km
    muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
    R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
    R_SOI_PL = 0.929*1e6; % km
elseif sim.ID_FLYBY == 4
    RPlanet_flyby = astroConstants(24); % Radius_Mars, km
    muPlanet_flyby = astroConstants(14); % mumarte, km^3/s^2
    R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
    R_SOI_PL = 0.578*1e6; %km
end
MJDGA1_dim = MJDGA1*sim.TU/(3600*24);
v_arr_GA1_dim = v_arr_GA1.*sim.DU./sim.TU;
v_dep_GA1_dim = v_dep_GA1.*sim.DU./sim.TU;
[delta_v_p_GA1, ~] = flyby(RPlanet_flyby, muPlanet_flyby, R_lim_from_planet, ...
                  MJDGA1_dim, v_arr_GA1_dim, v_dep_GA1_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_GA1), 'Not found')
    penalty_dv_GA1 = 10; 
end
% 1st leg - GA1 -> Ast 1
M_start_1st_leg = output_GA1.m(end); %  - sim.M_pods
[output_1] = NL_interpolator_of( r_GA1 , rA1 , v_dep_GA1 , vA1 , N_rev1 , TOF1 , M_start_1st_leg ,sim.PS.Isp , sim );
%     if ~isnan(output_1.Thrust(1,1)) && abs(output_1.t(end) - TOF1) < tol_TOF % if is not nan -> it's a valid solution
if abs(max(output_1.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg1 = 10*abs(max(output_1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end
%   2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end); %  - sim.M_pods
[output_2] = NL_interpolator_of( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
%     if ~isnan(output_2.Thrust(1,1)) && abs(output_2.t(end) - TOF2) < tol_TOF  % if is not nan -> it's a valid solution
if abs(max(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = 10*abs(max(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

% SC 2
% 0th leg - Earth -> GA2
[output_GA2] = NL_interpolator_of( r_EA , r_GA2 , v_dep , v_arr_GA2 , N_revGA2 , TOFGA_SC2 , sim.M2 ,sim.PS.Isp , sim );
if abs(max(output_GA2.T_magn)) > sim.max_Available_Thrust
    penalty_T_legGA2 = 10*abs(max(output_GA2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_GA2.t(end) - TOFGA_SC2) > tol_TOF
    penalty_TOF_legGA2 = abs(output_GA2.t(end) - TOFGA_SC2);
end
% gravity assist
v_dep_GA2 = v_GA2 + v_inf_magnGA2*[cos(el_GA2out)*cos(az_GA2out); cos(el_GA2out)*sin(az_GA2out); sin(el_GA2out)];
if sim.ID_FLYBY == 3
    RPlanet_flyby = astroConstants(23); % Radius_Earth, km
    muPlanet_flyby = astroConstants(13); % muEarth, km^3/s^2
    R_lim_from_planet = 500; % km, for earth is ok to avoid atmosphere
    R_SOI_PL = 0.929*1e6; % km
elseif sim.ID_FLYBY == 4
    RPlanet_flyby = astroConstants(24); % Radius_Mars, km
    muPlanet_flyby = astroConstants(14); % mumarte, km^3/s^2
    R_lim_from_planet = 200; % km, for mars is ok to avoid atmosphere
    R_SOI_PL = 0.578*1e6; %km
end
MJDGA2_dim = MJDGA2*sim.TU/(3600*24);
v_arr_GA2_dim = v_arr_GA2.*sim.DU./sim.TU;
v_dep_GA2_dim = v_dep_GA2.*sim.DU./sim.TU;
[delta_v_p_GA2, ~] = flyby(RPlanet_flyby, muPlanet_flyby, R_lim_from_planet, ...
                  MJDGA2_dim, v_arr_GA2_dim, v_dep_GA2_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_GA2), 'Not found')
    penalty_dv_GA2 = 10; 
end
% a_th leg - GA2 -> Ast_a
M_start_ath_leg = output_GA2.m(end); %  - sim.M_pods
[output_a] = NL_interpolator_of( r_GA2 , rAa , v_dep_GA2 , vAa , N_reva , TOFa , M_start_ath_leg ,sim.PS.Isp , sim );
%     if ~isnan(output_1.Thrust(1,1)) && abs(output_1.t(end) - TOF1) < tol_TOF % if is not nan -> it's a valid solution
if abs(max(output_a.T_magn)) > sim.max_Available_Thrust
    penalty_T_lega = 10*abs(max(output_a.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end
% b_th leg - Ast_a -> Ast_b
M_start_b_th_leg = output_a.m(end); %  - sim.M_pods
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , M_start_b_th_leg ,sim.PS.Isp , sim );
if abs(max(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = 10*abs(max(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

mass_fract_SC1 = (output_1.m(1) - output_2.m(end))/output_1.m(1);
mass_fract_SC2 = (output_a.m(1) - output_b.m(end))/output_a.m(1);

%     if mass_fract_SC1 > 0 && mass_fract_SC1 < 1 
%         if mass_fract_SC2 > 0 && mass_fract_SC2 < 1 
MF = max(mass_fract_SC1,mass_fract_SC2); % + mean(mf_sc1,mf_sc2), cosi sono piu o meno uguali
% disp('success')

obj_fun = MF + penalty_MAX_DURATION + ...
    penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb + ...
    penalty_T_legGA1 + penalty_TOF_legGA1 + penalty_dv_GA1 + ...
    penalty_T_legGA2 + penalty_TOF_legGA2 + penalty_dv_GA2;

end

