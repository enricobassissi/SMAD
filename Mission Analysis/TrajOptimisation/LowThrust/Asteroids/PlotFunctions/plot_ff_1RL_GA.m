function [output, r_encounter,r_departure,R_coasting,v_encounter, sol] = plot_ff_1RL_GA(x,sim,data, sol)

%% setting the inputs
MJD01 = x(1); % departure time from earth
TOFGA = x(2);
MJDPGA = MJD01 + TOFGA; % passage time GA
TOF1 = x(3);
MJDP1 = MJDPGA + TOF1; % arrival on ast 1
TOC1 = x(21);
MJDD1 = MJDP1 + TOC1; % departure from ast 1
TOF2 = x(4);
MJDP2 = MJDD1 + TOF2; % arrival on ast 2
TOC2 = x(22);
MJDD2 = MJDP2 + TOC2; % departure from ast 2
TOF3 = x(5);
MJDP3 = MJDD2 + TOF3; % arrival on ast 3
TOC3 = x(23);
MJDD3 = MJDP3 + TOC3; % departure from ast 3
TOF4 = x(6);
MJDP4 = MJDD3 + TOF4; % arrival on ast 4


% N REV GA
N_revGA = x(7);
% N REV1
N_rev1  = x(8);
% N REV2
N_rev2  = x(9);
% N REV3
N_rev3  = x(10);
% N REV4
N_rev4  = x(11);


% C3 launcher
v_inf_magn = x(13);
az = x(14);
elev = x(15);

% GA stuff IN
v_inf_magn2 = x(16);
az2 = x(17);
el2 = x(18);
% GA stuff out
az3 = x(19);
el3 = x(20);

% chosing which asteroid to visit
IDP = x(12); %index of permutation, the column of the Permutation Matrix of the asteroids
asteroid_1 = data.PermutationMatrix(IDP,1);
asteroid_2 = data.PermutationMatrix(IDP,2);
asteroid_3 = data.PermutationMatrix(IDP,3);
asteroid_4 = data.PermutationMatrix(IDP,4);

%%  Computing position and velocity of the planets in that days
% Departure from Earth
MJD01_dim = MJD01*sim.TU/(3600*24);
[kep_EA,ksun] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,ksun);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;

% GA
MJDPGA_dim = MJDPGA*sim.TU/(3600*24);
[kep_GA,ksun] = uplanet(MJDPGA_dim, sim.ID_FLYBY);
[r_GA, v_GA] = sv_from_coe(kep_GA,ksun);
r_GA = r_GA/sim.DU;
v_GA = v_GA/sim.DU*sim.TU;

% Arrival at 1st ast
MJDP1_dim = MJDP1*sim.TU/(3600*24);
[kep_ast_1] = uNEO2(MJDP1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
r1 = r1/sim.DU;
v1 = v1/sim.DU*sim.TU;


% Departure from 1st ast
MJDD1_dim = MJDD1*sim.TU/(3600*24);
[kep_ast_1d] = uNEO2(MJDD1_dim,asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1d, v1d] = sv_from_coe(kep_ast_1d,ksun); % km, km/s
r1d = r1d/sim.DU;
v1d = v1d/sim.DU*sim.TU;

% Coasting 1 on 1st ast 
%     time_int = [MJDP1_dim*24*3600, MJDD1_dim*24*3600];
%     y0 = [r1*sim.DU v1*sim.DU/sim.TU]; %km, km/s; 
%     options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
%     [~,rc1] = ode113(@rates, time_int, y0,options,'sun');
%     R_coasting.ast1 = rc1(:,[1 2 3]);
% R_coasting.ast1 = coasting_asteroids(MJDP1_dim,MJDD1_dim,asteroid_1);

% Arrival at 2nd ast
MJDP2_dim = MJDP2*sim.TU/(3600*24);
[kep_ast_2] = uNEO2(MJDP2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[r2, v2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
r2 = r2/sim.DU;
v2 = v2/sim.DU*sim.TU;

% Departure from 2nd ast
MJDD2_dim = MJDD2*sim.TU/(3600*24);
[kep_ast_2d] = uNEO2(MJDD2_dim,asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
[r2d, v2d] = sv_from_coe(kep_ast_2d,ksun); % km, km/s
r2d = r2d/sim.DU;
v2d = v2d/sim.DU*sim.TU;

% Coasting 2 on 2nd ast 
%     time_int2 = [MJDP2_dim*24*3600, MJDD2_dim*24*3600];
%     y02 = [r2*sim.DU v2*sim.DU/sim.TU]; %km, km/s; 
%     options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
%     [~,rc2] = ode113(@rates, time_int2, y02,options,'sun');
%     R_coasting.ast2 = rc2(:,[1 2 3]);
%      R_coasting.ast2 = coasting_asteroids(MJDP2_dim,MJDD2_dim,asteroid_2);

% Arrival at 3rd ast
MJDP3_dim = MJDP3*sim.TU/(3600*24);
[kep_ast_3] = uNEO2(MJDP3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
[r3, v3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
r3 = r3/sim.DU;
v3 = v3/sim.DU*sim.TU;

% Departure from 3rd ast
MJDD3_dim = MJDD3*sim.TU/(3600*24);
[kep_ast_3d] = uNEO2(MJDD3_dim,asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
[r3d, v3d] = sv_from_coe(kep_ast_3d,ksun); % km, km/s
r3d = r3d/sim.DU;
v3d = v3d/sim.DU*sim.TU;

% Coasting 3 on 3rd ast 
%     time_int3 = [MJDP3_dim*24*3600, MJDD3_dim*24*3600];
%     y03 = [r3*sim.DU v3*sim.DU/sim.TU]; %km, km/s; 
%     options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
%     [~,rc3] = ode113(@rates, time_int3, y03,options,'sun');
%     R_coasting.ast3 = rc3(:,[1 2 3]);
% R_coasting.ast3 = coasting_asteroids(MJDP3_dim,MJDD3_dim,asteroid_3);

% Arrival at 4th ast
MJDP4_dim = MJDP4*sim.TU/(3600*24);
[kep_ast_4] = uNEO2(MJDP4_dim,asteroid_4,data); % [km,-,rad,rad,rad,wrapped rad]
[r4, v4] = sv_from_coe(kep_ast_4,ksun); % km, km/s
r4 = r4/sim.DU;
v4 = v4/sim.DU*sim.TU;

%% launcher
v_launcher = v_inf_magn*[cos(elev)*cos(az); cos(elev)*sin(az); sin(elev)];
v_dep = v_EA + v_launcher;  %if parabolic escape (v_launcher = 0)

% Gravity assist
v_arr_GA = v_GA + v_inf_magn2*[cos(el2)*cos(az2); cos(el2)*sin(az2); sin(el2)];

%% NLI
% 1st leg - Earth -> GA
[sol.output_GA] = NL_interpolator( r_EA , r_GA , v_dep , v_arr_GA , N_revGA , TOFGA , sim.M ,sim.PS.Isp , sim );

% 2nd leg - GA -> Ast 1
v_dep_GA = v_GA + v_inf_magn2*[cos(el3)*cos(az3); cos(el3)*sin(az3); sin(el3)];

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

v_arr_GA_dim = v_arr_GA.*sim.DU./sim.TU;
v_dep_GA_dim = v_dep_GA.*sim.DU./sim.TU;
[delta_v_p,sol.rp] = flyby(RPlanet_flyby, muPlanet_flyby,R_lim_from_planet, ...
              MJDPGA_dim, v_arr_GA_dim, v_dep_GA_dim, sim.ID_FLYBY,R_SOI_PL);

M_after_GA = sol.output_GA.m(end);
[sol.output_1] = NL_interpolator( r_GA , r1 , v_dep_GA , v1 , N_rev1 , TOF1 , M_after_GA ,sim.PS.Isp , sim );


% 3rd leg - Ast1 -> Ast2
M_start_2nd_leg = sol.output_1.m(end)- sim.M_pods; %  - sim.M_pods;
[sol.output_2] = NL_interpolator( r1d , r2 , v1d , v2 , N_rev2 , TOF2 , M_start_2nd_leg ,sim.PS.Isp , sim );


% 4th leg - Ast2 -> Ast3
M_start_3rd_leg = sol.output_2.m(end)- sim.M_pods; %  - sim.M_pods;
[sol.output_3] = NL_interpolator( r2d , r3 , v2d , v3 , N_rev3 , TOF3 , M_start_3rd_leg ,sim.PS.Isp , sim );


% 5th leg - Ast3 -> Ast4
M_start_4th_leg = sol.output_3.m(end)- sim.M_pods; %  
[sol.output_4] = NL_interpolator( r3d , r4 , v3d , v4 , N_rev4 , TOF4 , M_start_4th_leg ,sim.PS.Isp , sim );


sol.mass_fract = (sol.output_GA.m(1) - sol.output_4.m(end))/sol.output_GA.m(1);

%put one after the other, all the thrust profiles
T_append = [sol.output_GA.Thrust(:,1),sol.output_GA.Thrust(:,2),sol.output_GA.Thrust(:,3);
            sol.output_1.Thrust(:,1),sol.output_1.Thrust(:,2),sol.output_1.Thrust(:,3);
            sol.output_2.Thrust(:,1),sol.output_2.Thrust(:,2),sol.output_2.Thrust(:,3);
            sol.output_3.Thrust(:,1),sol.output_3.Thrust(:,2),sol.output_3.Thrust(:,3);
            sol.output_4.Thrust(:,1),sol.output_4.Thrust(:,2),sol.output_4.Thrust(:,3)]; 

        
T = sqrt(T_append(:,1).^2 + T_append(:,3).^2);


                                   
% %% Output encounter states
r_encounter.EA = r_EA;
r_encounter.GA = r_GA;

r_encounter.ast1 = r1;
r_departure.ast1 = r1d;

% arr = MJD01 + (sol.output_GA.t(end) + sol.output_1.t(end))*sim.TU/(3600*24);
% dep = arr + TOC1*sim.TU/(3600*24);
arr = MJD01 + (TOFGA + TOF1)*sim.TU/(3600*24);
dep = arr + TOC1*sim.TU/(3600*24);
time_int = [arr*24*3600, dep*24*3600];
y0 = [r1*sim.DU v1*sim.DU/sim.TU]; %km, km/s; 
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[~,rc1] = ode113(@rates, time_int, y0,options,'sun');
AU = astroConstants(2);
R_coasting.ast1 = rc1(:,[1 2 3])./AU;
% R_coasting.ast1 = coasting_asteroids(MJDP1_dim,MJDD1_dim,asteroid_1);
%R_coasting.ast1 = coasting_asteroids(arr,dep,asteroid_1);

r_encounter.ast2 = r2;
r_departure.ast2 = r2d;


arr = MJD01 + (sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end))*sim.TU/(3600*24);
dep = arr + TOC2*sim.TU/(3600*24);
time_int = [arr*24*3600, dep*24*3600];
y0 = [r2*sim.DU v2*sim.DU/sim.TU]; %km, km/s; 
[~,rc2] = ode113(@rates, time_int, y0,options,'sun');
%R_coasting.ast2 = rc2(:,[1 2 3])./AU;
R_coasting.ast2 = coasting_asteroids(arr,dep,asteroid_2);


r_encounter.ast3 = r3;
r_departure.ast3 = r3d;

arr = MJD01 + (sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end) + TOC2 + sol.output_1.t(end))*sim.TU/(3600*24);
dep = arr + TOC3*sim.TU/(3600*24);
time_int = [arr*24*3600, dep*24*3600];
y0 = [r3*sim.DU v3*sim.DU/sim.TU]; %km, km/s; 
[~,rc3] = ode113(@rates, time_int, y0,options,'sun');
R_coasting.ast3 = rc3(:,[1 2 3])./AU;

r_encounter.ast4 = r4;

v_encounter.EA = v_EA;
v_encounter.GA = v_GA;
v_encounter.ast1 = v1;
v_departure.ast1 = v1;
v_encounter.ast2 = v2;
v_encounter.ast3 = v3;
v_encounter.ast4 = v4;
% add the others if you need


%% Output

output.t            = [sol.output_GA.t;
                       sol.output_GA.t(end) + sol.output_1.t; 
                       sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t; 
                       sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end) + TOC2 + sol.output_3.t; 
                       sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end) + TOC2 + sol.output_3.t(end) + TOC3 + sol.output_4.t]; 
                   
output.m            = [sol.output_GA.m; sol.output_1.m; sol.output_2.m; sol.output_3.m; sol.output_4.m] ;
output.Thrust       = [sol.output_GA.Thrust; sol.output_1.Thrust; sol.output_2.Thrust; sol.output_3.Thrust; sol.output_4.Thrust];
output.a            = [sol.output_GA.a; sol.output_1.a; sol.output_2.a; sol.output_3.a; sol.output_4.a];

output.r.GA         = sol.output_GA.r;
output.r.leg1       = sol.output_1.r; 
output.r.leg2       = sol.output_2.r;
output.r.leg3       = sol.output_3.r;
output.r.leg4       = sol.output_4.r;
output.theta.GA     = sol.output_GA.theta; 
output.theta.leg1   = sol.output_1.theta; 
output.theta.leg2   = sol.output_2.theta;
output.theta.leg3   = sol.output_3.theta;
output.theta.leg4   = sol.output_4.theta;
output.z.GA         = sol.output_GA.z;
output.z.leg1       = sol.output_1.z;
output.z.leg2       = sol.output_2.z;
output.z.leg3       = sol.output_3.z;
output.z.leg4       = sol.output_4.z;
output.Href.GA      = sol.output_GA.Href;
output.Href.leg1    = sol.output_1.Href;
output.Href.leg2    = sol.output_2.Href;
output.Href.leg3    = sol.output_3.Href;
output.Href.leg4    = sol.output_4.Href;

output.tEnd.LegGA2 = sol.output_GA.t(end)*sim.TU/86400; 
output.tEnd.Leg12 = (sol.output_GA.t(end) + sol.output_1.t(end))*sim.TU/86400;
output.tEnd.Leg22 = (sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end))*sim.TU/86400;
output.tEnd.Leg32 = (sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end) + TOC2)*sim.TU/86400;
output.tEnd.Leg42 = (sol.output_GA.t(end) + sol.output_1.t(end) + TOC1 + sol.output_2.t(end) + TOC2 + sol.output_3.t(end) + TOC3)*sim.TU/86400;

end

