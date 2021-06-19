function [output, r_encounter, v_encounter, sol] = plot_ff_2RL_GA_ball_simplified(x,sim,data,sol,power_propulsion_data)
%{ 

input variable vector
x = [
% Departure dates (1)
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
% -- ID Permutation (13)]
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

% -- N REV1
N_rev1 = round(x(9));
% N REV2
N_rev2 = round(x(12));

% OUT
az_out_GA1 = x(6);
el_out_GA1 = x(7);

% -- DSM SC1
rD_mag1 = x(3);
thetaD1 = x(4);


%% choosing which asteroid to visit
% 1ST SPACECRAFT ASTEROID OBJECTIVES
% IDP1 = x(29); %index of permutation, the column of the Permutation Matrix of the asteroids
% asteroid_1 = data.PermutationMatrix(IDP1,1);
% asteroid_2 = data.PermutationMatrix(IDP1,2);
asteroid_1 = sol.asteroid_1;
asteroid_2 = sol.asteroid_2;

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

%% FUEL CONSUMPTIONS %%
tol_TOF = 1; % 1 TU means approx 60 days
penalty_dv_GA1 = 0;
penalty_T_leg1 = 0; penalty_T_leg2 = 0;
penalty_TOF_leg1 = 0; penalty_TOF_leg2 = 0; 

Isp_secondary_prop = 255; % ADN
g0 = 9.81;

%% Both SC togheter
% ballistic stuff before GA
% DV calculation with lambert
% Earth -> DSM
[~,~,~,~,V_EA_DSM_dep,V_EA_DSM_arr,~,~] = lambertMR(r_EA,rDSM1,TOF_DSM1,sim.mu,0,0,0,0);
V_EA_DSM_dep = V_EA_DSM_dep';
V_EA_DSM_arr = V_EA_DSM_arr';
sqrt_C3_requested = sqrt((V_EA_DSM_dep(1)-v_EA(1))^2+(V_EA_DSM_dep(2)-v_EA(2))^2+(V_EA_DSM_dep(3)-v_EA(3))^2);
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
v_inf_minus = V_DSM_GA_arr - v_GA1;
v_inf_plus = norm(v_inf_minus)*[cos(el_out_GA1)*cos(az_out_GA1); cos(el_out_GA1)*sin(az_out_GA1); sin(el_out_GA1)];
v_dep_GA1 = v_inf_plus + v_GA1; 
Mass_at_beginning_LT = sim.M1*exp(-dV_DSM_GA*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));
[output_1] = NL_interpolator_of( r_GA1 , rA1 , v_dep_GA1 , vA1 , N_rev1 , TOF1 , Mass_at_beginning_LT ,sim.PS.Isp , sim );
if max(abs(output_1.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg1 = max(abs(output_1.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_1.t(end) - TOF1) > tol_TOF
    penalty_TOF_leg1 = abs(output_1.t(end) - TOF1);
end

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if TOF_DSM1+TOFGA_1+TOF1+CT1+TOF2 > max_duration
    penalty_MAX_DURATION = TOF_DSM1+TOFGA_1+TOF1+CT1+TOF2 - max_duration; % 12 years max mission time 
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

% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end) - sim.M_pods; % yes pods here, released on Ast1
[output_2] = NL_interpolator_of( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
if max(abs(output_2.T_magn)) > sim.max_Available_Thrust
    penalty_T_leg2 = max(abs(output_2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_2.t(end) - TOF2) > tol_TOF
    penalty_TOF_leg2 = abs(output_2.t(end) - TOF2);
end

%% masses and obj function
Isp = sim.PS.Isp*sim.TU;

% SC1
mass_depleted_DSM1 = sim.M1 - Mass_at_beginning_LT;
mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);

sol.mass_depleted_DSM1 = mass_depleted_DSM1;
sol.mass_depleted_Leg1 = mass_depleted_Leg1;
sol.mass_depleted_Leg2 = mass_depleted_Leg2;

sol.tot_mass_depleted_SC1 = mass_depleted_DSM1+mass_depleted_Leg1+mass_depleted_Leg2;
sol.mass_dry_and_pods_SC1 = sim.M1 - sol.tot_mass_depleted_SC1;
mass_fract_SC1 = sol.tot_mass_depleted_SC1/sim.M1;

sol.dV_DSM1 = dV_DSM_GA*sim.DU/sim.TU*1e3;
sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

% gravity assist quantities
sol.delta_V_p1 = delta_v_p_1;
sol.dV_gained_GA1 = (sqrt((v_dep_GA1(1)-V_DSM_GA_arr(1))^2+(v_dep_GA1(2)-V_DSM_GA_arr(2))^2+(v_dep_GA1(3)-V_DSM_GA_arr(3))^2) - sol.delta_V_p1)*sim.DU/sim.TU;

penalty_MF_unfeasible = 0;
if mass_fract_SC1 <= 0 || mass_fract_SC1 >= 1 
    penalty_MF_unfeasible = 100;
end

sol.mass_fract_SC1 = mass_fract_SC1; 

sol.obj_fun = mass_fract_SC1 + penalty_MAX_DURATION + penalty_MF_unfeasible + ...
    20*dV_DSM_GA + penalty_dV_C3 + ...
    penalty_dv_GA1 + ...
    20*(penalty_T_leg1 + penalty_T_leg2) + ...
    penalty_TOF_leg1 + penalty_TOF_leg2 + ....
    pen_dV_DSM;
    
%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.DSM1 = rDSM1;
r_encounter.GA1 = r_GA1;
r_encounter.astA1 = rA1;
r_encounter.astD1 = rD1;
r_encounter.astA2 = rA2;

v_encounter.EA = v_EA;
v_encounter.dep_EA_DSM1 = V_EA_DSM_dep;
v_encounter.dep_DSM1_GA1 = V_DSM_GA_dep;
v_encounter.GA1 = v_GA1;
v_encounter.astA1 = vA1;
v_encounter.astD1 = vD1;
v_encounter.astA2 = vA2;

%% Porcherie
% t_span_EA_DSM1 = linspace(0,TOF_DSM1*sim.TU,sim.n_sol);
% t_span_DSM1_GA1 = linspace(TOF_DSM1*sim.TU,(TOF_DSM1+TOFGA_1)*sim.TU,sim.n_sol);
% t_span_CT1 = linspace((TOF_DSM1+TOFGA_1)*sim.TU+output_1.t(end),(TOF_DSM1+TOFGA_1)*sim.TU+output_1.t(end)+CT1,sim.n_sol);
% t_span_EA_DSM2 = linspace(0,TOF_DSM2*sim.TU,sim.n_sol);
% t_span_DSM2_GA2 = linspace(TOF_DSM2*sim.TU,(TOF_DSM2+TOFGA_a)*sim.TU,sim.n_sol);
% t_span_CTa = linspace((TOF_DSM2+TOFGA_a)*sim.TU+output_a.t(end),(TOF_DSM2+TOFGA_a)*sim.TU+output_a.t(end)+CTa,sim.n_sol);
% only LT plotting
t_span_CT1 = linspace(output_1.t(end),output_1.t(end)+CT1,sim.n_sol);
mCT1 = ones(sim.n_sol,1).*output_1.m(end);
TCT1 = zeros(sim.n_sol,3); 

output.t_SC1            = [output_1.t; t_span_CT1'; t_span_CT1(end)+output_2.t];
output.m_SC1            = [output_1.m; mCT1; output_2.m];
output.Thrust_SC1       = [output_1.Thrust; TCT1; output_2.Thrust];
output.T_magn_SC1       = sqrt(output.Thrust_SC1(:,1).^2 + output.Thrust_SC1(:,3).^2);
output.a_SC1            = [output_1.a; output_2.a];
output.r.leg1       = output_1.r; 
output.r.leg2       = output_2.r;
output.theta.leg1   = output_1.theta; 
output.theta.leg2   = output_2.theta;
output.z.leg1       = output_1.z;
output.z.leg2       = output_2.z;
output.Href.leg1    = output_1.Href;
output.Href.leg2    = output_2.Href;

output.t1           = output_1.t;
output.CT1          = t_span_CT1';
output.t2           = output_2.t;

%% --
sol.T_1 = output_1.Thrust;
sol.T_2 = output_2.Thrust;

end