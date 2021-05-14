function [output, r_encounter, v_encounter] = plot_ff_ea_2ast_LT_soo_NLI(x,sim,data)
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

% N REV1
N_rev1 = x(4);
% N REV2
N_rev2 = x(5);
% N REV3
N_rev3 = x(13);

% C3 launcher
v_inf_magn = x(7);
alpha = x(8);
beta = x(9);

% chosing which asteroid to visit
IDP = x(6); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
% asteroid_4 = data.PermutationMatrix(IDP,4);

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

%% launcher
v_launcher = v_inf_magn*[cos(alpha)*cos(beta); sin(alpha)*cos(beta); sin(beta)] ;
v_dep = v_EA + v_launcher;  %since parabolic escape (vinf = 0)

%% 1st leg - Earth -> Ast 1
[output_1] = NL_interpolator( r_EA , rA1 , v_dep , vA1 , N_rev1 , TOF1 , sim.M ,sim.PS.Isp , sim );

% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end);
[output_2] = NL_interpolator( rD1 , rA2 , vD1 , vA2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );

% 3rd leg - Ast2 -> Ast3
M_start_3rd_leg = output_2.m(end);
[output_3] = NL_interpolator( rD2 , rA3 , vD2 , vA3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );

%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.astA1 = rA1;
r_encounter.astD1 = rD1;
r_encounter.astA2 = rA2;
r_encounter.astD2 = rD2;
r_encounter.astA3 = rA3;

v_encounter.EA = v_EA;
v_encounter.astA1 = vA1;
v_encounter.astD1 = vD1;
v_encounter.astA2 = vA2;
v_encounter.astD2 = vD2;
v_encounter.astA3 = vA3;

%% Porcherie
t_span_CT1 = linspace(output_1.t(end),output_1.t(end)+CT1,sim.n_sol);
t_span_CT2 = linspace(t_span_CT1(end)+output_2.t(end),t_span_CT1(end)+output_2.t(end)+CT2,sim.n_sol);
mCT1 = ones(sim.n_sol,1).*output_1.m(end);
mCT2 = ones(sim.n_sol,1).*output_2.m(end);
TCT1 = ones(sim.n_sol,3).*output_1.Thrust(end);
TCT2 = ones(sim.n_sol,3).*output_2.Thrust(end);
%% Output
output.t            = [output_1.t; t_span_CT1'; t_span_CT1(end)+output_2.t; t_span_CT2'; ...
                       t_span_CT2(end) + output_3.t];
output.m            = [output_1.m; mCT1; output_2.m; mCT2; output_3.m] ;
output.Thrust       = [output_1.Thrust; TCT1; output_2.Thrust; TCT2; output_3.Thrust];
output.a            = [output_1.a; output_2.a; output_3.a];
output.r.leg1       = output_1.r; 
output.r.leg2       = output_2.r;
output.r.leg3       = output_3.r;
output.theta.leg1   = output_1.theta; 
output.theta.leg2   = output_2.theta;
output.theta.leg3   = output_3.theta;
output.z.leg1       = output_1.z;
output.z.leg2       = output_2.z;
output.z.leg3       = output_3.z;
output.Href.leg1    = output_1.Href;
output.Href.leg2    = output_2.Href;
output.Href.leg3    = output_3.Href;

end

