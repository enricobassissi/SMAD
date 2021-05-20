function [output, r_encounter, v_encounter] = plot_ff_1ast_to_another_ast_LT_soo_NLI(x,sim,data)
% setting the input times
MJDP1 = x(1);
TOF = x(2);
MJDP2 = MJDP1 + TOF;

% N REV
N_rev = x(3);

% chosing which asteroid to visit
IDP = x(4); %index of permutation, the column of the Permutation Matrix of the asteroids
IDP2 = x(5);
% asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_1 = data.asteroid_names(IDP);
TF = contains(data.asteroid_names,asteroid_1);
data_elements_matrix_2nd_ast = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2nd_ast, ~] = ...
            sequences_local_pruning(data_elements_matrix_2nd_ast, 1);
asteroid_2 = PermutationMatrix_2nd_ast(IDP2);
% asteroid_2 = data.PermutationMatrix(IDP,2);
% asteroid_3 = data.PermutationMatrix(IDP,3);
% asteroid_4 = data.PermutationMatrix(IDP,4);

% Computing position and velocity of the planets in that days
% % Departure from Earth
% MJD01_dim = MJD01*sim.TU/(3600*24);
% [kep_EA,ksun] = uplanet(MJD01_dim, 3);
% [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
% r_EA = r_EA/sim.DU;
% v_EA = v_EA/sim.DU*sim.TU;

% passage at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,sim.mu_dim); % km, km/s
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;

% passage at 2nd ast
MJDP2_dim = MJDP2*sim.TU/(3600*24);
[kep_ast_2] = uNEO2(MJDP2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[r2, v2] = sv_from_coe(kep_ast_2,sim.mu_dim); % km, km/s
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% v_launcher = v_inf_magn*[cos(el)*cos(az); cos(el)*sin(az); sin(el)] ;
% v_dep = v1 + v_launcher;  %since parabolic escape (vinf = 0)
% varr = v2; %since we want to rendezvous

[output] = NL_interpolator( r1 , r2 , v1 , v2 , N_rev , TOF , sim.M ,sim.PS.Isp , sim );

r_encounter.ast1 = r1;
r_encounter.ast2 = r2;

v_encounter.ast1 = v1;
v_encounter.ast2 = v2;

end

