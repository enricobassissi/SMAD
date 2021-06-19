function [output, r_encounter, sol] = plot_ff_2FL(x,sim,data, sol)
%% setting the input times
MJD01 = x(1); % departure time from earth
TOF1 = x(2);
MJDP1 = MJD01 + TOF1; % passage time on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2; % passage ast 2

TOFa = x(4);
MJDPa = MJD01 + TOFa; % passage ast 3
TOFb = x(5);
MJDPb = MJDPa + TOFb; % passage ast 4

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
el = x(14);

% chosing which asteroids to visit

% 1st SPACECRAFT ASTEROIDS OBJECTIVES
IDP = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);

% 2nd SPACECRAFT ASTEROIDS OBJECTIVES
IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
asteroid_a = PermutationMatrix_2SC(IDP2,1);
asteroid_b = PermutationMatrix_2SC(IDP2,2);


% asteroid1 flyby
v_inf_ast1_magn = x(15);
az_ast1 = x(16);
el_ast1 = x(17);

% asteroid2 flyby
v_inf_ast2_magn = x(18);
az_ast2 = x(19);
el_ast2 = x(20);

% asteroida flyby
v_inf_asta_magn = x(21);
az_asta = x(22);
el_asta = x(23);

% asteroidb flyby
v_inf_astb_magn = x(24);
az_astb = x(25);
el_astb = x(26);

%% Computing position and velocity of planets/asteroids in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;


% PASSAGE at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% passage at 2nd ast
MJDP2_dim = MJDP2*sim.TU/(3600*24);
[kep_ast_2] = uNEO2(MJDP2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% passage at a ast
MJDPa_dim = MJDPa*sim.TU/(3600*24);
[kep_ast_a] = uNEO2(MJDPa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[ra, va] = sv_from_coe(kep_ast_a,ksun); % km, km/s
ra = ra/sim.DU;
va = va/sim.DU*sim.TU;

% passage at b ast
MJDPb_dim = MJDPb*sim.TU/(3600*24);
[kep_ast_b] = uNEO2(MJDPb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rb, vb] = sv_from_coe(kep_ast_b,ksun); % km, km/s
rb = rb/sim.DU;
vb = vb/sim.DU*sim.TU;


%% Relative velocity -> heliocentric velocity
% Launcher departure variable
v_extralaunch = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)];
v_dep = v_EA + v_extralaunch;  %if parabolic escape (v_extralaunch = 0)

% Flyby ast 1
v_rel_ast1 = v_inf_ast1_magn* [cos(el_ast1)*cos(az_ast1); cos(el_ast1)*sin(az_ast1); sin(el_ast1)];
v_abs_ast1 = v1 + v_rel_ast1;

% Flyby ast 2
v_rel_ast2 = v_inf_ast2_magn* [cos(el_ast2)*cos(az_ast2); cos(el_ast2)*sin(az_ast2); sin(el_ast2)];
v_abs_ast2 = v2 + v_rel_ast2;

% Flyby ast a
v_rel_asta = v_inf_asta_magn* [cos(el_asta)*cos(az_asta); cos(el_asta)*sin(az_asta); sin(el_asta)];
v_abs_asta = va + v_rel_asta;

% Flyby ast b
v_rel_astb = v_inf_astb_magn* [cos(el_astb)*cos(az_astb); cos(el_astb)*sin(az_astb); sin(el_astb)];
v_abs_astb = vb + v_rel_astb;

%% NLI
Isp = sim.PS.Isp*sim.TU;
g0 = sim.g0/sim.TU^2*(1000*sim.DU);
% SC1
% 1st leg - Earth -> ast1
[output_1] = NL_interpolator( r_EA , r1 , v_dep , v_abs_ast1, N_rev1 , TOF1 , sim.M1 ,sim.PS.Isp , sim );

sol.mass_depleted_Leg1 = output_1.m(1) - output_1.m(end);
sol.dV_associated_Leg1 = -g0*Isp*log(output_1.m(end)/output_1.m(1)); % -ve*ln(m_final/m_initial)

% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end) - sim.M_pods; %  
[output_2] = NL_interpolator( r1 , r2 , v_abs_ast1 , v_abs_ast2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );

sol.mass_depleted_Leg2 = output_2.m(1) - output_2.m(end);
sol.dV_associated_Leg2 = -g0*Isp*log(output_2.m(end)/output_2.m(1)); % -ve*ln(m_final/m_initial)

% sol.mass_fract_SC1 = (output_1.m(1) - output_2.m(end))/output_1.m(1);
% there are the pods! and they are mass ok but not propellant
sol.tot_mass_depleted_SC1 = sol.mass_depleted_Leg1+sol.mass_depleted_Leg2;
sol.mass_dry_and_pods_SC1 = output_1.m(1) - sol.tot_mass_depleted_SC1;
sol.mass_fract_SC1 = (output_1.m(1) - sol.mass_dry_and_pods_SC1)/output_1.m(1);

sol.T_1 = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3)];
sol.T_2 = [output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3)];

T_append1 = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3);
             output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3)];

T1 = sqrt(T_append1(:,1).^2 + T_append1(:,3).^2);


% SC2
% 1st leg - Earth -> asta
[output_a] = NL_interpolator( r_EA , ra , v_dep , v_abs_asta, N_reva , TOFa , sim.M2 ,sim.PS.Isp , sim );

sol.mass_depleted_Lega = output_a.m(1) - output_a.m(end);
sol.dV_associated_Lega = -g0*Isp*log(output_a.m(end)/output_a.m(1)); % -ve*ln(m_final/m_initial)

% 2nd leg - Asta -> Astb
M_start_b_leg = output_a.m(end) - sim.M_pods; % 
[output_b] = NL_interpolator( ra , rb , v_abs_asta , v_abs_astb , N_revb , TOFb , M_start_b_leg ,sim.PS.Isp , sim );

sol.mass_depleted_Legb = output_b.m(1) - output_b.m(end);
sol.dV_associated_Legb = -g0*Isp*log(output_b.m(end)/output_b.m(1)); % -ve*ln(m_final/m_initial)

% sol.mass_fract_SC2 = (output_a.m(1) - output_b.m(end))/output_a.m(1);
% there are the pods! and they are mass ok but not propellant
sol.tot_mass_depleted_SC2 = sol.mass_depleted_Lega+sol.mass_depleted_Legb;
sol.mass_dry_and_pods_SC2 = output_a.m(1) - sol.tot_mass_depleted_SC2;
sol.mass_fract_SC2 = (output_a.m(1) - sol.mass_dry_and_pods_SC2)/output_a.m(1);

sol.T_a = [output_a.Thrust(:,1),output_a.Thrust(:,2),output_a.Thrust(:,3)];
sol.T_b = [output_b.Thrust(:,1),output_b.Thrust(:,2),output_b.Thrust(:,3)];

T_append2 = [output_a.Thrust(:,1),output_a.Thrust(:,2),output_a.Thrust(:,3);
             output_b.Thrust(:,1),output_b.Thrust(:,2),output_b.Thrust(:,3)];

T2 = sqrt(T_append2(:,1).^2 + T_append2(:,3).^2);


%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.ast1 = r1;
r_encounter.ast2 = r2;
r_encounter.asta = ra;
r_encounter.astb = rb;

%% Output
output.t_SC1            = [output_1.t; output_1.t(end)+output_2.t];
output.t_SC2            = [output_a.t; output_a.t(end)+output_b.t];

output.m_SC1            = [output_1.m; output_2.m];
output.m_SC2            = [output_a.m; output_b.m];

output.Thrust_SC1       = T_append1;
output.T_magn_SC1       = T1;
output.Thrust_SC2       = T_append2;
output.T_magn_SC2       = T2;

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
output.t2           = output_2.t;
output.ta           = output_a.t;
output.tb           = output_b.t;

end

