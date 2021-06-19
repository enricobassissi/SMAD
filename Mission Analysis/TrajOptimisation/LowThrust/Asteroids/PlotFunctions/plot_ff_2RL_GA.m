function [output, r_encounter, v_encounter, sol] = plot_ff_2RL_GA(x,sim,data, sol)
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
N REV GA (18)
vinf (GRAVITY ASSIST) (19)
azimuth in (PRE-GRAVITY ASSIST) (20)
elevation (PRE-GRAVITY ASSIST) (21)
azimuth out (POST-GRAVITY ASSIST)(22)
elevation (POST-GRAVITY ASSIST)(23)

% GA Stuff SC2
TOFGA (24)
N REV GA (25)
vinf (GRAVITY ASSIST) (26)
azimuth in (PRE-GRAVITY ASSIST) (27)
elevation (PRE-GRAVITY ASSIST) (28)
azimuth out (POST-GRAVITY ASSIST)(29)
elevation (POST-GRAVITY ASSIST)(30)]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

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
TOFGA_a = x(24);
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
% N REV GA1
N_rev_GA1 = x(18);
% N REV GAa
N_rev_GAa = x(25);

% C3 launcher
v_inf_magn = x(12);
az = x(13);
elev = x(14);

% GA stuff SC1
% IN
v_inf_magn_GA1 = x(19);
az_in_GA1 = x(20);
el_in_GA1 = x(21);
% OUT
az_out_GA1 = x(22);
el_out_GA1 = x(23);

% GA stuff SC2
% IN
v_inf_magn_GAa = x(26);
az_in_GAa = x(27);
el_in_GAa = x(28);
% OUT
az_out_GAa = x(29);
el_out_GAa = x(30);

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
v_GA1 = v_GA1/sim.DU*sim.TU;

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

% GA SC1
MJDPGAa_dim = MJDPGA_a*sim.TU/(3600*24);
[kep_GAa,ksun] = uplanet(MJDPGAa_dim, sim.ID_FLYBY);
[r_GAa, v_GAa] = sv_from_coe(kep_GAa,ksun);
r_GAa = r_GAa/sim.DU;
v_GAa = v_GAa/sim.DU*sim.TU;

% passage at b_th ast
MJDAb_dim = MJDAb*sim.TU/(3600*24);
[kep_ast_Ab] = uNEO2(MJDAb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,ksun); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;

%% Launcher departure variable
v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

% Gravity assist SC1
v_arr_GA1 = v_GA1 + v_inf_magn_GA1*[cos(el_in_GA1)*cos(az_in_GA1); cos(el_in_GA1)*sin(az_in_GA1); sin(el_in_GA1)];

% Gravity assist SC2
v_arr_GAa = v_GAa + v_inf_magn_GAa*[cos(el_in_GAa)*cos(az_in_GAa); cos(el_in_GAa)*sin(az_in_GAa); sin(el_in_GAa)];

%% NLI %%
tol_TOF = 0.5; % 1 TU means approx 60 days
penalty_T_legGA1 = 0; penalty_TOF_legGA1 = 0; penalty_dv_GA1 = 0;
penalty_T_legGAa = 0; penalty_TOF_legGAa = 0; penalty_dv_GAa = 0;
penalty_T_leg1 = 0; penalty_T_leg2 = 0; penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; penalty_TOF_lega = 0; penalty_TOF_legb = 0;

%% SC1 
% 1st leg - Earth -> Ast 1
[output_1] = NL_interpolator_of( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , sim.M1 ,sim.PS.Isp , sim );
%     if ~isnan(output_1.Thrust(1,1)) && abs(output_1.t(end) - TOF1) < tol_TOF % if is not nan -> it's a valid solution
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    %penalty_T_leg1 = abs(max(output_1.T_magn)) - sim.max_Available_Thrust;
    penalty_T_leg1 = max(abs(output_1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end

% 2nd leg - Ast1 -> GA SC1
M_start_GA1_leg = output_1.m(end) - sim.M_pods; % yes pods here, released on Ast1
[output_GA1] = NL_interpolator_of( rD1 , r_GA1 , vD1 , v_arr_GA1 , N_rev_GA1 , TOFGA_1 , M_start_GA1_leg ,sim.PS.Isp , sim );
%     if ~isnan(output_2.Thrust(1,1)) && abs(output_2.t(end) - TOF2) < tol_TOF  % if is not nan -> it's a valid solution
if max(abs(output_GA1.T_magn)) > sim.max_Available_Thrust
    %penalty_T_leg2 = abs(max(output_2.T_magn)) - sim.max_Available_Thrust;
    penalty_T_legGA1 = max(abs(output_GA1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_GA1.t(end) - TOFGA_1) > tol_TOF
    penalty_TOF_legGA1 = abs(output_GA1.t(end) - TOFGA_1);
end

v_dep_GA1 = v_GA1 + v_inf_magn_GA1*[cos(el_out_GA1)*cos(az_out_GA1); cos(el_out_GA1)*sin(az_out_GA1); sin(el_out_GA1)];
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
    penalty_dv_GA1 = 20;
end

% 3rd leg - GA SC1 -> Ast 2
M_after_GA1 = output_GA1.m(end); % no pods here because it's the GA
[output_2] = NL_interpolator_of( r_GA1 , rA2 , v_dep_GA1 , vA2 , N_rev2 , TOF2 , M_after_GA1 ,sim.PS.Isp , sim );
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = abs(max(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

%% SC 2
% a_th leg - Earth -> Ast_a
[output_a] = NL_interpolator_of( r_EA , rAa , v_dep , vAa , N_reva , TOFa , sim.M2 ,sim.PS.Isp , sim );
%     if ~isnan(output_a.Thrust(1,1)) && abs(output_a.t(end) - TOFa) < tol_TOF % if is not nan -> it's a valid solution
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
    %penalty_T_lega = 10*abs(max(output_a.T_magn)) - sim.max_Available_Thrust;
    penalty_T_lega = 10*max(abs(output_a.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_a.t(end) - TOFa) > tol_TOF
    penalty_TOF_lega = abs(output_a.t(end) - TOFa);
end

% GAa Leg - Ast_a -> GAa
M_start_GAa_leg = output_a.m(end) - sim.M_pods; % yes pods here because you release mass on Ast_a
[output_GAa] = NL_interpolator_of( rDa , r_GAa , vDa , v_arr_GAa , N_rev_GAa , TOFGA_a , M_start_GAa_leg ,sim.PS.Isp , sim );
if max(abs(output_GAa.T_magn)) > sim.max_Available_Thrust
    %penalty_T_leg2 = abs(max(output_2.T_magn)) - sim.max_Available_Thrust;
    penalty_T_legGAa = max(abs(output_GAa.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_GAa.t(end) - TOFGA_a) > tol_TOF
    penalty_TOF_legGAa = abs(output_GAa.t(end) - TOFGA_a);
end

v_dep_GAa = v_GAa + v_inf_magn_GAa*[cos(el_out_GAa)*cos(az_out_GAa); cos(el_out_GAa)*sin(az_out_GAa); sin(el_out_GAa)];
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
    penalty_dv_GAa = 20;
end

% b_th Leg - GAa -> Ast_b
M_start_b_th_leg = output_GAa.m(end); % no pods here because it's the GA
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , M_start_b_th_leg ,sim.PS.Isp , sim );
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    %penalty_T_legb = 10*abs(max(output_b.T_magn)) - sim.max_Available_Thrust;
    penalty_T_legb = 10*max(abs(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

%% masses and obj function
Isp = sim.PS.Isp*sim.TU;
g0 = sim.g0/sim.TU^2*(1000*sim.DU);

% SC1
sol.mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
sol.mass_depleted_LegGA1 = output_GA1.m(1) - output_GA1.m(end);
sol.mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);

sol.tot_mass_depleted_SC1 = sol.mass_depleted_Leg1+sol.mass_depleted_LegGA1+sol.mass_depleted_Leg2;
sol.mass_dry_and_pods_SC1 = output_1.m(1) - sol.tot_mass_depleted_SC1;
sol.mass_fract_SC1 = (output_1.m(1) - sol.mass_dry_and_pods_SC1)/output_1.m(1);

sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_LegGA1 = -g0*Isp*log(output_GA1.m(end)/output_GA1.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

% SC2
sol.mass_depleted_Lega = output_a.m(1) - output_a.m(end);
sol.mass_depleted_LegGAa = output_GAa.m(1) - output_GAa.m(end);
sol.mass_depleted_Legb = output_b.m(1) - output_b.m(end);

sol.tot_mass_depleted_SC2 = sol.mass_depleted_Lega+sol.mass_depleted_LegGAa+sol.mass_depleted_Legb;
sol.mass_dry_and_pods_SC2 = output_a.m(1) - tot_mass_depleted_SC2;
sol.mass_fract_SC2 = (output_a.m(1) - sol.mass_dry_and_pods_SC2)/output_a.m(1);

sol.dV_associated_Lega = -g0*Isp*log(output_a.m(end)/output_a.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_LegGAa = -g0*Isp*log(output_GAa.m(end)/output_GAa.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_Legb = -g0*Isp*log(output_b.m(end)/output_b.m(1)); % -ve*ln(m_final/m_initial)

avg_mass_fraction = (sol.mass_fract_SC1+sol.mass_fract_SC2)/2;
MF = max(sol.mass_fract_SC1,sol.mass_fract_SC2) + abs(sol.mass_fract_SC1 - avg_mass_fraction) + ...
        abs(sol.mass_fract_SC2 - avg_mass_fraction); % cosi sono piu o meno uguali

sol.obj_fun = MF + penalty_MAX_DURATION + penalty_dv_GA1 + penalty_dv_GAa + ...
    20*(penalty_T_leg1 + penalty_T_leg2 + penalty_T_lega + penalty_T_legb + ...
    penalty_T_legGA1 + penalty_T_legGAa) + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + penalty_TOF_lega + penalty_TOF_legb + ...
    penalty_TOF_legGA1 + penalty_TOF_legGAa;


%% Output encounter states
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
TCT1 = ones(sim.n_sol,3).*output_1.Thrust(end); %% oppure è thrust nulla???
TCTa = ones(sim.n_sol,3).*output_a.Thrust(end);

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

