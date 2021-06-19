function [output, r_encounter, v_encounter, sol] = plot_ff_1RL(x,sim,data, sol)
%% setting the input times
MJD01 = x(1); % departure time from earth
TOF1 = x(2);
MJDA1 = MJD01 + TOF1; % arrival time on ast 1
CT1 = x(10); % coasting time
MJDD1 = MJDA1 + CT1; % departure from ast 1
TOF2 = x(3);
MJDA2 = MJDD1 + TOF2; % arrival ast 2
CT2 = x(11); % coasting time 2
MJDD2 = MJDA2 + CT2; % departure on ast 2
TOF3 = x(12);
MJDA3 = MJDD2 + TOF3; % arrival ast 3
CT3 = x(14); % coasting time 3
MJDD3 = MJDA3 + CT3; % departure on ast 3
TOF4 = x(15);
MJDA4 = MJDD3 + TOF4; % arrival ast 4

% N REV1
N_rev1 = x(4);
% N REV2
N_rev2 = x(5);
% N REV3
N_rev3 = x(13);
% N REV4
N_rev4 = x(16);

% C3 launcher
v_inf_magn = x(7);
az = x(8);
el = x(9);

% chosing which asteroid to visit
IDP = x(6); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);

%% Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% arrival at 1st ast
MJDA1_dim = MJDA1*sim.TU/(3600*24);
[kep_ast_A1] = uNEO2(MJDA1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rA1, vA1] = sv_from_coe(kep_ast_A1,ksun); % km, km/s
rA1 = rA1/sim.DU;
vA1 = vA1/sim.DU*sim.TU;

% departure from 1st ast
MJDD1_dim = MJDD1*sim.TU/(3600*24);
[kep_ast_D1] = uNEO2(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[rD1, vD1] = sv_from_coe(kep_ast_D1,ksun); % km, km/s
rD1 = rD1/sim.DU;
vD1 = vD1/sim.DU*sim.TU;

% arrival at 2nd ast
MJDA2_dim = MJDA2*sim.TU/(3600*24);
[kep_ast_A2] = uNEO2(MJDA2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rA2, vA2] = sv_from_coe(kep_ast_A2,ksun); % km, km/s
rA2 = rA2/sim.DU;
vA2 = vA2/sim.DU*sim.TU;

% departure from 2nd ast
MJDD2_dim = MJDD2*sim.TU/(3600*24);
[kep_ast_D2] = uNEO2(MJDD2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[rD2, vD2] = sv_from_coe(kep_ast_D2,ksun); % km, km/s
rD2 = rD2/sim.DU;
vD2 = vD2/sim.DU*sim.TU;

% arrival at 3rd ast
MJDA3_dim = MJDA3*sim.TU/(3600*24);
[kep_ast_A3] = uNEO2(MJDA3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
[rA3, vA3] = sv_from_coe(kep_ast_A3,ksun); % km, km/s
rA3 = rA3/sim.DU;
vA3 = vA3/sim.DU*sim.TU;

% departure from 3rd ast
MJDD3_dim = MJDD3*sim.TU/(3600*24);
[kep_ast_D3] = uNEO2(MJDD3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
[rD3, vD3] = sv_from_coe(kep_ast_D3,ksun); % km, km/s
rD3 = rD3/sim.DU;
vD3 = vD3/sim.DU*sim.TU;

% arrival at 4th ast
MJDA4_dim = MJDA4*sim.TU/(3600*24);
[kep_ast_A4] = uNEO2(MJDA4_dim,asteroid_4,data); % [km,-,rad,rad,rad,wrapped rad]
[rA4, vA4] = sv_from_coe(kep_ast_A4,ksun); % km, km/s
rA4 = rA4/sim.DU;
vA4 = vA4/sim.DU*sim.TU;

%% launcher
v_launcher = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)
sol.v_launcher = v_launcher;

%% NLI
Isp = sim.PS.Isp*sim.TU;
g0 = sim.g0/sim.TU^2*(1000*sim.DU);

% 1st leg - Earth -> Ast 1
[output_1] = NL_interpolator( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , sim.M ,sim.PS.Isp , sim );

sol.mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)

% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end) - sim.M_pods;
[output_2] = NL_interpolator( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );

sol.mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

% 3rd leg - Ast2 -> Ast3
M_start_3rd_leg = output_2.m(end) - sim.M_pods;
[output_3] = NL_interpolator( rD2 , rA3 , vD2 , vA3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );

sol.mass_depleted_Leg3 = output_3.m(1) - output_3.m(end);
sol.dV_associated_Leg3 = -g0*Isp*log(output_3.m(end)/output_3.m(1)); % -ve*ln(m_final/m_initial)

% 4th leg - Ast3 -> Ast4
M_start_4th_leg = output_3.m(end) - sim.M_pods;
[output_4] = NL_interpolator( rD3 , rA4 , vD3 , vA4 , N_rev4 , TOF4 , M_start_4th_leg ,sim.PS.Isp , sim );

sol.mass_depleted_Leg4 = output_4.m(1) - output_4.m(end);
sol.dV_associated_Leg4 = -g0*Isp*log(output_4.m(end)/output_4.m(1)); % -ve*ln(m_final/m_initial)

%% Extract quantities from output struct
% sol.mass_fract = (output_1.m(1) - output_4.m(end))/output_1.m(1); % now
% there are the pods! and they are mass ok but not propellant
sol.tot_mass_depleted = sol.mass_depleted_Leg1+sol.mass_depleted_Leg2+sol.mass_depleted_Leg3+sol.mass_depleted_Leg4;
sol.mass_dry_and_pods = output_1.m(1) - sol.tot_mass_depleted;
sol.mass_fract = (output_1.m(1) - sol.mass_dry_and_pods)/output_1.m(1);

sol.T_1 = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3)];
sol.T_2 = [output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3)];
sol.T_3 = [output_3.Thrust(:,1),output_3.Thrust(:,2),output_3.Thrust(:,3)];
sol.T_4 = [output_4.Thrust(:,1),output_4.Thrust(:,2),output_4.Thrust(:,3)];

%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.astA1 = rA1;
r_encounter.astD1 = rD1;
r_encounter.astA2 = rA2;
r_encounter.astD2 = rD2;
r_encounter.astA3 = rA3;
r_encounter.astD3 = rD3;
r_encounter.astA4 = rA4;

v_encounter.EA = v_EA;
v_encounter.astA1 = vA1;
v_encounter.astD1 = vD1;
v_encounter.astA2 = vA2;
v_encounter.astD2 = vD2;
v_encounter.astA3 = vA3;
v_encounter.astD3 = vD3;
v_encounter.astA4 = vA4;

%% Porcherie
t_span_CT1 = linspace(output_1.t(end),output_1.t(end)+CT1,sim.n_sol);
t_span_CT2 = linspace(t_span_CT1(end)+output_2.t(end),t_span_CT1(end)+output_2.t(end)+CT2,sim.n_sol);
t_span_CT3 = linspace(t_span_CT2(end)+output_3.t(end),t_span_CT2(end)+output_3.t(end)+CT3,sim.n_sol);
mCT1 = ones(sim.n_sol,1).*output_1.m(end);
mCT2 = ones(sim.n_sol,1).*output_2.m(end);
mCT3 = ones(sim.n_sol,1).*output_3.m(end);
TCT1 = ones(sim.n_sol,3).*output_1.Thrust(end);
TCT2 = ones(sim.n_sol,3).*output_2.Thrust(end);
TCT3 = ones(sim.n_sol,3).*output_3.Thrust(end);

%% Output
output.t            = [output_1.t; t_span_CT1'; t_span_CT1(end)+output_2.t; t_span_CT2'; ...
                       t_span_CT2(end) + output_3.t ; t_span_CT3'; t_span_CT3(end) + output_4.t];
output.m            = [output_1.m; mCT1; output_2.m; mCT2; output_3.m; mCT3; output_4.m] ;
output.Thrust       = [output_1.Thrust; TCT1; output_2.Thrust; TCT2; output_3.Thrust; TCT3; output_4.Thrust];
output.a            = [output_1.a; output_2.a; output_3.a; output_4.a];
output.r.leg1       = output_1.r; 
output.r.leg2       = output_2.r;
output.r.leg3       = output_3.r;
output.r.leg4       = output_4.r;
output.theta.leg1   = output_1.theta; 
output.theta.leg2   = output_2.theta;
output.theta.leg3   = output_3.theta;
output.theta.leg4   = output_4.theta;
output.z.leg1       = output_1.z;
output.z.leg2       = output_2.z;
output.z.leg3       = output_3.z;
output.z.leg4       = output_4.z;
output.Href.leg1    = output_1.Href;
output.Href.leg2    = output_2.Href;
output.Href.leg3    = output_3.Href;
output.Href.leg4    = output_4.Href;

%% single legs outputs
output.t1 = output_1.t;
output.t2 = output_2.t;
output.t3 = output_3.t;
output.t4 = output_4.t;

output.tEnd.Leg11 = output.t1(end)*sim.TU/86400; %  A1
output.tEnd.Leg12 = output.tEnd.Leg11+CT1*sim.TU/86400; %  D1
output.tEnd.Leg21 = output.tEnd.Leg12+output_2.t(end)*sim.TU/86400; %  A2
output.tEnd.Leg22 = output.tEnd.Leg21+CT2*sim.TU/86400; % ct2 on ast 2
output.tEnd.Leg31 = output.tEnd.Leg22+output_3.t(end)*sim.TU/86400; % A3
output.tEnd.Leg32 = output.tEnd.Leg31+CT3*sim.TU/86400; %  D3
output.tEnd.Leg41 = output.tEnd.Leg32+output_4.t(end)*sim.TU/86400; %  A4
end

