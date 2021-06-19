function [output, r_encounter, v_encounter, sol] = plot_ff_2RL_GA_ball_simplified_JP(x,sim,data,sol)
%{ 

input variable vector
% x = [
% rP1 (1)
% thetaP1 (2)
% TOF_P1 (3), from earth dep to P1 time
% rP2 (4)
% thetaP2 (5)
% TOF_P2 (6), from P1 to P2
% rDSM (7)
% thetaDSM (8)
% TOF_DSM (9), from P2 to DSM
% rP3 (10)
% thetaP3 (11)
% TOF_P3 (12), from DSM1 to P3
% rP4 (13)
% thetaP4 (14)
% TOF_P4 (15), from P3 to P4
% TOF_GA (16)
% azimuth out (POST-GRAVITY ASSIST)(17)
% elevation (POST-GRAVITY ASSIST)(18)
% TOFa (19)
% N REV a (20)
% Coasting time a (21)
% TOFb (22)
% N REV b (23)
% % -- ID Permutation (24)]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)
MJD01 = x(25);
TOF_P1 = x(3); % from earth dep to P1 time
TOF_P2 = x(6); % from P1 to P2
TOF_DSM = x(9); % from P2 to DSM
TOF_P3 = x(12); % from DSM1 to P3
TOF_P4 = x(15); % from P3 to P4
TOF_GAa = x(16);
TOFa = x(19);
CTa = x(21);
TOFb = x(22);

rP1_magn = x(1);
thetaP1 = x(2);

rP2_magn = x(4);
rel_thetaP2 = x(5);

rDSM_magn = x(7);
rel_thetaDSM = x(8);

rP3_magn = x(10);
rel_thetaP3 = x(11);

rP4_magn = x(13);
rel_thetaP4 = x(14);

az_out_GAa = x(17);
el_out_GAa = x(18);

N_reva = x(20);
N_revb = x(23);

iP1 = x(26);
iP2 = x(27);
iDSM = x(28);
iP3 = x(29);
iP4 = x(30);

%% choosing which asteroid to visit
IDP = x(24);
asteroid_a = data.PermutationMatrix(IDP,1);
asteroid_b = data.PermutationMatrix(IDP,2);

% arbitrarly fixed, they're nice in sequence
% asteroid_a = "2010UJ";
% asteroid_b = "2016CK137";

%% FUEL CONSUMPTIONS %%
tol_TOF = 1; % 1 TU means approx 60 days

penalty_T_P2 = 0; penalty_TOF_P2 = 0;
penalty_T_P4 = 0; penalty_TOF_P4 = 0;
penalty_dv_GAa = 0;
penalty_T_lega = 0; penalty_T_legb = 0; 
penalty_TOF_lega = 0; penalty_TOF_legb = 0;

Isp_secondary_prop = 255; % ADN
g0 = 9.81;

%% ------ Times and Positions
MJDP1 = MJD01 + TOF_P1; % passage time P1
MJDP2 = MJDP1 + TOF_P2; % passage time P2
MJD_DSM = MJDP2 + TOF_DSM; % passage time DSM
MJDP3 = MJD_DSM + TOF_P3; % passage time P3
MJDP4 = MJDP3 + TOF_P4; % passage time P4
MJDPGA_a = MJDP4 + TOF_GAa; % passage at GA
MJDAa = MJDPGA_a + TOFa; % mjd2000 arrive of SC2 on ast a
MJDDa = MJDAa + CTa; % departure from ath asteroid
MJDAb = MJDDa + TOFb; % mjd2000 passage of SC2 on ast b

% Points Passages
rP1 = rP1_magn*[cos(thetaP1)*cos(iP1); sin(thetaP1)*cos(iP1); sin(iP1)];
thetaP2 = thetaP1 + rel_thetaP2;
rP2 = rP2_magn*[cos(thetaP2)*cos(iP2); sin(thetaP2)*cos(iP2); sin(iP2)];
thetaDSM = thetaP2 + rel_thetaDSM;
rDSM = rDSM_magn*[cos(thetaDSM)*cos(iDSM); sin(thetaDSM)*cos(iDSM); sin(iDSM)];
thetaP3 = thetaDSM + rel_thetaP3;
rP3 = rP3_magn*[cos(thetaP3)*cos(iP3); sin(thetaP3)*cos(iP3); sin(iP3)];
thetaP4 = thetaP3 + rel_thetaP4;
rP4 = rP4_magn*[cos(thetaP4)*cos(iP4); sin(thetaP4)*cos(iP4); sin(iP4)];

% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

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

%% ------------ Manoeuvres --------------- %
% Earth Dep -> P1, Lambert
[~,~,~,~,V_EA_P1_dep,V_EA_P1_arr,~,~] = lambertMR(r_EA,rP1,TOF_P1,sim.mu,0,0,0,0);
V_EA_P1_dep = V_EA_P1_dep';
V_EA_P1_arr = V_EA_P1_arr';
dV_EA_P1 = sqrt((V_EA_P1_dep(1)-v_EA(1))^2+(V_EA_P1_dep(2)-v_EA(2))^2+(V_EA_P1_dep(3)-v_EA(3))^2);
if dV_EA_P1 < sim.dV_DSM_max % vinf that the launcher can give max 
    pen_dV_EA_P1 = 0;
else
    pen_dV_EA_P1 = 1e3*abs(dV_EA_P1 - sim.dV_DSM_max); % penalty like, but not discard a priori
end

% P2 -> DSM, Lambert
[~,~,~,~,V_P2_DSM_dep,V_P2_DSM_arr,~,~] = lambertMR(rP2,rDSM,TOF_DSM,sim.mu,0,0,0,0);
V_P2_DSM_dep = V_P2_DSM_dep';
V_P2_DSM_arr = V_P2_DSM_arr';

% P1 -> P2, LowThrust
% NREV = 0, i want it to be short
M_start_P1_P2_leg = sim.M1*exp(-dV_EA_P1*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));
% M_start_P1_P2_leg = sim.M2;
[output_P2] = NL_interpolator_of( rP1 , rP2 , V_EA_P1_arr , V_P2_DSM_dep , 0 , TOF_P2 , M_start_P1_P2_leg ,sim.PS.Isp , sim );
% [ output_P2 ] = CW_LowLambert_of( rP1 ,rP2 , V_EA_P1_arr , V_P2_DSM_dep , 0 , TOF_P2 ,M_start_P1_P2_leg ,3 , 3 , sim.PS ,sim )

if max(abs(output_P2.T_magn)) > sim.max_Available_Thrust
    penalty_T_P2 = max(abs(output_P2.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_P2.t(end) - TOF_P2) > tol_TOF
    penalty_TOF_P2 = abs(output_P2.t(end) - TOF_P2);
end

% DSM -> P3, Lambert
[~,~,~,~,V_DSM_P3_dep,V_DSM_P3_arr,~,~] = lambertMR(rDSM,rP3,TOF_P3,sim.mu,0,0,0,0);
V_DSM_P3_dep = V_DSM_P3_dep';
V_DSM_P3_arr = V_DSM_P3_arr';
dV_DSM_P3 = sqrt((V_DSM_P3_dep(1)-V_P2_DSM_arr(1))^2+(V_DSM_P3_dep(2)-V_P2_DSM_arr(2))^2+(V_DSM_P3_dep(3)-V_P2_DSM_arr(3))^2);
if dV_DSM_P3 < sim.dV_DSM_max % vinf that the launcher can give max 
    pen_dV_DSM_P3 = 0;
else
    pen_dV_DSM_P3 = 1e3*abs(dV_DSM_P3 - sim.dV_DSM_max); % penalty like, but not discard a priori
end

% P4 -> GA
[~,~,~,~,V_P4_GA_dep,V_P4_GA_arr,~,~] = lambertMR(rP4,r_GAa,TOF_GAa,sim.mu,0,0,0,0);
V_P4_GA_dep = V_P4_GA_dep';
V_P4_GA_arr = V_P4_GA_arr';

% P3 -> P4, LowThrust
% NREV = 0, i want it to be short
M_start_P3_P4_leg = output_P2.m(end)*exp(-dV_DSM_P3*(sim.DU/sim.TU)*1e3/(g0*Isp_secondary_prop));
% M_start_P3_P4_leg = output_P2.m(end);
[output_P4] = NL_interpolator_of( rP3 , rP4 , V_DSM_P3_arr , V_P4_GA_dep , 0 , TOF_P4 , M_start_P3_P4_leg ,sim.PS.Isp , sim );
if max(abs(output_P4.T_magn)) > sim.max_Available_Thrust
    penalty_T_P4 = max(abs(output_P4.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_P4.t(end) - TOF_P4) > tol_TOF
    penalty_TOF_P4 = abs(output_P4.t(end) - TOF_P4);
end

% ath leg - GA -> Ast a
v_inf_minus = V_P4_GA_arr - v_GAa;
v_inf_plus = norm(v_inf_minus)*[cos(el_out_GAa)*cos(az_out_GAa); cos(el_out_GAa)*sin(az_out_GAa); sin(el_out_GAa)];
v_dep_GAa = v_inf_plus + v_GAa; 
Mass_at_beginning_LT2 = output_P4.m(end);
[output_a] = NL_interpolator_of( r_GAa , rAa , v_dep_GAa , vAa , N_reva , TOFa , Mass_at_beginning_LT2 ,sim.PS.Isp , sim );
if max(abs(output_a.T_magn)) > sim.max_Available_Thrust
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
v_arr_GAa_dim = V_P4_GA_arr.*sim.DU./sim.TU;
v_dep_GAa_dim = v_dep_GAa.*sim.DU./sim.TU;
[delta_v_p_a,~] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
                  MJDPGAa_dim, v_arr_GAa_dim, v_dep_GAa_dim, sim.ID_FLYBY, R_SOI_PL);
if strcmp(string(delta_v_p_a), 'Not found')
    penalty_dv_GAa = 10;
end

% bth leg - Asta -> Astb
M_start_bth_leg = output_a.m(end) - sim.M_pods; % yes pods here, released on Ast1
[output_b] = NL_interpolator_of( rDa , rAb , vDa , vAb , N_revb , TOFb , M_start_bth_leg ,sim.PS.Isp , sim );
if max(abs(output_b.T_magn)) > sim.max_Available_Thrust
    penalty_T_legb = max(abs(output_b.T_magn)) - sim.max_Available_Thrust;
end
if abs(output_b.t(end) - TOFb) > tol_TOF
    penalty_TOF_legb = abs(output_b.t(end) - TOFb);
end

%% masses and obj function
Isp = sim.PS.Isp*sim.TU;

sol.mass_depleted_EA_P1 = sim.M1 - M_start_P1_P2_leg;
sol.mass_depleted_P1_P2 = output_P2.m(1) - output_P2.m(end);
sol.mass_depleted_DSM = output_P2.m(end) - M_start_P3_P4_leg;
sol.mass_depleted_P3_P4 = output_P4.m(1) - output_P4.m(end);
sol.mass_depleted_Lega = output_a.m(1) - output_a.m(end);
sol.mass_depleted_Legb = output_b.m(1) - output_b.m(end);

sol.tot_mass_depleted_SC2 = sol.mass_depleted_EA_P1+sol.mass_depleted_P1_P2+sol.mass_depleted_DSM+...
    sol.mass_depleted_P3_P4+sol.mass_depleted_Lega+sol.mass_depleted_Legb;
sol.mass_fract_SC2 = sol.tot_mass_depleted_SC2/sim.M1;

sol.dV_DSM2 = dV_DSM_P3*sim.DU/sim.TU*1e3;
sol.dV_associated_Lega = -g0*Isp*log(output_a.m(end)/output_a.m(1)); % -ve*ln(m_final/m_initial)
sol.dV_associated_Legb = -g0*Isp*log(output_b.m(end)/output_b.m(1)); % -ve*ln(m_final/m_initial)

% gravity assist quantities
sol.delta_V_pa = delta_v_p_a;
sol.dV_gained_GAa = (sqrt((v_dep_GAa(1)-V_P4_GA_arr(1))^2+(v_dep_GAa(2)-V_P4_GA_arr(2))^2+(v_dep_GAa(3)-V_P4_GA_arr(3))^2) - sol.delta_V_pa)*sim.DU/sim.TU;

penalty_MF_unfeasible = 0;
if sol.mass_fract_SC2 <= 0 || sol.mass_fract_SC2 >= 1
    penalty_MF_unfeasible = 100;
end

max_duration = 12*365*(3600*24)/sim.TU;
penalty_MAX_DURATION = 0;
if TOF_P1+TOF_P2+TOF_DSM+TOF_P3+TOF_P4+TOF_GAa+TOFa+CTa+TOFb > max_duration
    penalty_MAX_DURATION = TOF_P1+TOF_P2+TOF_DSM+TOF_P3+TOF_P4+TOF_GAa+TOFa+CTa+TOFb - max_duration; % 12 years max mission time 
end

sol.obj_fun = sol.mass_fract_SC2 + penalty_MAX_DURATION + penalty_MF_unfeasible + ...
    20*(pen_dV_EA_P1 + pen_dV_DSM_P3) + penalty_dv_GAa + ...
    20*(penalty_T_P2 + penalty_T_P4 + penalty_T_lega + penalty_T_legb) + ...
    penalty_TOF_P2 + penalty_TOF_P4 + penalty_TOF_lega + penalty_TOF_legb;
    
%% Output encounter states
r_encounter.DSM2 = rDSM;
r_encounter.GAa = r_GAa;
r_encounter.astAa = rAa;
r_encounter.astDa = rDa;
r_encounter.astAb = rAb;

v_encounter.dep_DSM2_GAa = V_DSM_P3_dep;
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
t_span_CTa = linspace(output_a.t(end),output_a.t(end)+CTa,sim.n_sol);
mCTa = ones(sim.n_sol,1).*output_a.m(end);
TCTa = zeros(sim.n_sol,3); 

output.t_SC2            = [output_a.t; t_span_CTa'; t_span_CTa(end)+output_b.t];
output.m_SC2            = [output_a.m; mCTa; output_b.m];
output.Thrust_SC2       = [output_a.Thrust; TCTa; output_b.Thrust];
output.T_magn_SC2       = sqrt(output.Thrust_SC2(:,1).^2 + output.Thrust_SC2(:,3).^2);
output.a_SC2            = [output_a.a; output_b.a];
output.r.lega       = output_a.r; 
output.r.legb       = output_b.r;
output.theta.lega   = output_a.theta; 
output.theta.legb   = output_b.theta;
output.z.lega       = output_a.z;
output.z.legb       = output_b.z;
output.Href.lega    = output_a.Href;
output.Href.legb    = output_b.Href;

output.ta           = output_a.t;
output.CTa          = t_span_CTa';
output.tb           = output_b.t;

%% --
sol.T_a = output_a.Thrust;
sol.T_b = output_b.Thrust;

sol.T_max = max(output.T_magn_SC2);

end