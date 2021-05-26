function [output, r_encounter, sol] = plot_ff_1FL_Ale(x,sim,data, sol)
%% setting the input times
MJD01 = x(1); % departure time from earth
TOF1 = x(2);
MJDP1 = MJD01 + TOF1; % passage time on ast 1
TOF2 = x(3);
MJDP2 = MJDP1 + TOF2; % passage ast 2
TOF3 = x(4);
MJDP3 = MJDP2 + TOF3; % passage ast 3
TOF4 = x(5);
MJDP4 = MJDP3 + TOF4; % passage ast 4


% N REV1
N_rev1 = x(6);
% N REV2
N_rev2 = x(7);
% N REV3
N_rev3 = x(8);
% N REV4
N_rev4 = x(9);

% C3 launcher
v_inf_magn = x(11);
az = x(12);
el = x(13);

% GA stuff IN
v_inf_magn2 = x(15);
az2 = x(16);
el2 = x(17);
% GA stuff out
az3 = x(18);
el3 = x(19);

% chosing which asteroid to visit
IDP = x(10); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);


% asteroid1 flyby
v_inf_ast1_magn = x(14);
az_ast1 = x(15);
el_ast1 = x(16);

% asteroid2 flyby
v_inf_ast2_magn = x(17);
az_ast2 = x(18);
el_ast2 = x(19);

% asteroid3 flyby
v_inf_ast3_magn = x(20);
az_ast3 = x(21);
el_ast3 = x(22);

% asteroid4 flyby
v_inf_ast4_magn = x(23);
az_ast4 = x(24);
el_ast4 = x(25);

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

% passage at 3rd ast
MJDP3_dim = MJDP3*sim.TU/(3600*24);
[kep_ast_3] = uNEO2(MJDP3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
[r3, v3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
r3 = r3/sim.DU;
v3 = v3/sim.DU*sim.TU;

% passage at 4th ast
MJDP4_dim = MJDP4*sim.TU/(3600*24);
[kep_ast_4] = uNEO2(MJDP4_dim,asteroid_4,data); % [km,-,rad,rad,rad,wrapped rad]
[r4, v4] = sv_from_coe(kep_ast_4,ksun); % km, km/s
r4 = r4/sim.DU;
v4 = v4/sim.DU*sim.TU;

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

% Flyby ast 3
v_rel_ast3 = v_inf_ast3_magn* [cos(el_ast3)*cos(az_ast3); cos(el_ast3)*sin(az_ast3); sin(el_ast3)];
v_abs_ast3 = v3 + v_rel_ast3;

% Flyby ast 4
v_rel_ast4 = v_inf_ast4_magn* [cos(el_ast4)*cos(az_ast4); cos(el_ast4)*sin(az_ast4); sin(el_ast4)];
v_abs_ast4 = v4 + v_rel_ast4;

%% NLI

% 1ST leg - EA -> Ast 1
[output_1] = NL_interpolator( r_EA , r1 , v_dep , v_abs_ast1, N_rev1 , TOF1 , sim.M ,sim.PS.Isp , sim );

% 3rd leg - Ast1 -> Ast2
M_start_2nd_leg = output_1.m(end); %  - sim.M_pods;
[output_2] = NL_interpolator( r1 , r2 , v_abs_ast1 , v_abs_ast2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );

% 4th leg - Ast2 -> Ast3
M_start_3rd_leg = output_2.m(end); %  - sim.M_pods;
[output_3] = NL_interpolator( r2 , r3 , v_abs_ast2 , v_abs_ast3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );

% 5th leg - Ast3 -> Ast4
M_start_4th_leg = output_3.m(end); %  - sim.M_pods;
[output_4] = NL_interpolator( r3 , r4 , v_abs_ast3 , v_abs_ast4 , N_rev4 , TOF4 , M_start_4th_leg ,sim.PS.Isp , sim );

sol.mass_fract = (output_1.m(1) - output_4.m(end))/output_1.m(1);


T_append = [output_1.Thrust(:,1),output_1.Thrust(:,2),output_1.Thrust(:,3);
            output_2.Thrust(:,1),output_2.Thrust(:,2),output_2.Thrust(:,3);
            output_3.Thrust(:,1),output_3.Thrust(:,2),output_3.Thrust(:,3);
            output_4.Thrust(:,1),output_4.Thrust(:,2),output_4.Thrust(:,3)];

sol.T = sqrt(T_append(:,1).^2 + T_append(:,3).^2);



%% Output encounter states
r_encounter.EA = r_EA;
r_encounter.ast1 = r1;
r_encounter.ast2 = r2;
r_encounter.ast3 = r3;
r_encounter.ast4 = r4;


%% Output
output.t            = [output_1.t; output_1.t(end)+output_2.t; ...
                       output_1.t(end)+output_2.t(end)+output_3.t; ...
                       output_1.t(end)+output_2.t(end)+output_3.t(end)+output_4.t];
                   
output.m            = [output_1.m; output_2.m; output_3.m; output_4.m] ;

output.Thrust       = [output_1.Thrust; output_2.Thrust; ...
                       output_3.Thrust; output_4.Thrust];
                   
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


output.t1           = output_1.t;
output.t2           = output_2.t;
output.t3           = output_3.t;
output.t4           = output_4.t;

end

