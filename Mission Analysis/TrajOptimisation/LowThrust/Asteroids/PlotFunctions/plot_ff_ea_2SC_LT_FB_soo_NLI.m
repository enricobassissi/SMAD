function [output, r_encounter, v_encounter, sol] = plot_ff_ea_2SC_LT_FB_soo_NLI(x,sim,data, sol)
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
     (14) elevation]
%}

% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% setting the input times
MJD01 = x(1); % departure time for both the sc

% 1st spacecraft characteristic times
TOF1 = x(2); % tof sc1 to 1st asteroid
MJDP1 = MJD01 + TOF1; % mjd2000 passage of 1st sc on ast 1
TOF2 = x(3); % tof sc1 to 2nd asteroid
MJDP2 = MJDP1 + TOF2; % mjd2000 passage of 1st sc on ast 2

% 2nd spacecraft characteristic times
TOFa = x(4); % tof sc2 to 1st asteroid
MJDPa = MJD01 + TOFa; 
TOFb =  x(5); % tof sc2 to 2nd asteroid
MJDPb = MJDPa + TOFb; 

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

%% choosing which asteroid to visit
% 1ST SPACECRAFT ASTEROID OBJECTIVES
IDP1 = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP1,1);
asteroid_2 = data.PermutationMatrix(IDP1,2);

% 2ND SPACECRAFT ASTEROID OBJECTIVES
IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
asteroid_sequence = [asteroid_1,asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_element_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning(data_elements_matrix_2SC, data.p_number);
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

% passage at a_th ast
MJDPa_dim = MJDPa*sim.TU/(3600*24);
[kep_ast_a] = uNEO2(MJDPa_dim,asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[ra, va] = sv_from_coe(kep_ast_a,ksun); % km, km/s
ra = ra/sim.DU;
va = va/sim.DU*sim.TU;

% passage at b_th ast
MJDPb_dim = MJDPb*sim.TU/(3600*24);
[kep_ast_b] = uNEO2(MJDPb_dim,asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rb, vb] = sv_from_coe(kep_ast_b,ksun); % km, km/s
rb = rb/sim.DU;
vb = vb/sim.DU*sim.TU;

%% Launcher departure variable
%     z = r .* sin(elev);
%     rcoselev = r .* cos(elev);
%     x = rcoselev .* cos(az);
%     y = rcoselev .* sin(az);
v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_extra = 0)

%% 1st leg - Earth -> Ast 1
% SC1
% 1st leg - Earth -> Ast1
[output_1] = NL_interpolator( r_EA , r1 , v_dep , v1 , N_rev1 , TOF1 , sim.M1 ,sim.PS.Isp , sim );
sol.output_1 = output_1;
% 2nd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end); %  - sim.M_pods
[output_2] = NL_interpolator( r1 , r2 , v1 , v2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );
sol.output_2 = output_2;
% SC 2
[output_a] = NL_interpolator( r_EA , ra , v_dep , va , N_reva , TOFa , sim.M2 ,sim.PS.Isp , sim );
sol.output_a = output_a;
% a_th leg - Earth -> Ast_a
M_start_b_th_leg = output_a.m(end); %  - sim.M_pods
[output_b] = NL_interpolator( ra , rb , va , vb , N_revb , TOFb , M_start_b_th_leg ,sim.PS.Isp , sim );
sol.output_b = output_b;
% mass fractions
output.mass_fract_SC1 = (output_1.m(1) - output_2.m(end))/output_1.m(1);
output.mass_fract_SC2 = (output_a.m(1) - output_b.m(end))/output_a.m(1);

%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.ast1 = r1;
r_encounter.ast2 = r2;
r_encounter.asta = ra;
r_encounter.astb = rb;

v_encounter.EA = v_EA;
v_encounter.ast1 = v1;
v_encounter.ast2 = v2;
v_encounter.asta = va;
v_encounter.astb = vb;

%% Output
output.t_SC1            = [output_1.t; output_1.t(end)+output_2.t];
output.t_SC2            = [output_a.t; output_a.t(end)+output_b.t];
output.m_SC1            = [output_1.m; output_2.m];
output.m_SC2            = [output_a.m; output_b.m];
output.Thrust_SC1       = [output_1.Thrust; output_2.Thrust];
output.T_magn_SC1       = sqrt(output.Thrust_SC1(:,1).^2 + output.Thrust_SC1(:,3).^2);
output.Thrust_SC2       = [output_a.Thrust; output_b.Thrust];
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
output.t2           = output_2.t;
output.ta           = output_a.t;
output.tb           = output_b.t;

end

