%% --------------------------------------------------------------------- %%
%% --------------- Earth Ast1 Ast2 Ast3 Ast4 Transfer ------------------ %%
%% ---------------------- ARCH 1+4, LT GA FB --------------------------- %%
%% --------------------- FREEDOM AT EACH AST --------------------------- %%
%% ------------------------------ SOO ---------------------------------- %%
%% --------------------------------------------------------------------- %%
%% Setup for default options
set(0, 'DefaultTextFontSize', 20)
set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultLegendFontSize', 20)
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.8)
format short

%% Initializing the Environment
clear; close all; clc;

% Palette ESA
colors = [0    50   71;... % (1) DEEP SPACE
          207  29   57;... % (2) EXCITE RED +1
          0    142  122;... % (3) PURE TEAL +1
          251  171  24;... % (4) ENLIGHT YELLOW
          244  121  32;... % (5) ENLIGHT YELLOW +1
          150  1    54;... % (6) EXCITE RED +2
          167  85   52;... % (7) ENLIGHT YELLOW +2
          0    97   158;... % (8) TRUSTY AZURE +1
          30   51   120;... % (9) TRUSTY AZURE +2
          0    103  98;... % (10) PURE TEAL +2
          51   94   111;... % (11) DEEP SPACE -1
          0    0    0]./255; % (12) BLACK

sim.case_name = 'ARCH ID 6: LOW THRUST RV GA ON EVERY ASTEROID';

%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%% Asteroids
% load('data_elements_matrix_9.mat')
% data.p_number = 4;
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix, data.p_number);
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

%load('data_processed_9_4.mat')

load('data_processed_42_475.mat')

%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 200;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
sim.M = 100; % SC wet mass [kg]
sim.M_pods = 5; % mass of the payloads+landing stuff [kg]
sim.ID_FLYBY = 4; % flyby planet

sim.c = 1e+3; %case 12: 1e+3
sim.max_Available_Thrust = 0.01; %case 12: 0.02

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld = [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOFGA (14)
bound.TOFGA_min = 0.5*365*3600*24/sim.TU; %%
bound.TOFGA_max = 3*365*3600*24/sim.TU; 
% TOF1 (2)
bound.TOF1_min = 0.3*365*3600*24/sim.TU; 
bound.TOF1_max = 3*365*3600*24/sim.TU; 
% TOF2 (3)
bound.TOF2_min = 0.3*365*3600*24/sim.TU; 
bound.TOF2_max = 3*365*3600*24/sim.TU; 
% TOF3 (4)
bound.TOF3_min = 0.3*365*3600*24/sim.TU; 
bound.TOF3_max = 3*365*3600*24/sim.TU; 
% TOF4 (5)
bound.TOF4_min = 0.3*365*3600*24/sim.TU; 
bound.TOF4_max = 3*365*3600*24/sim.TU; 
% N REV 1 (6)
bound.N_REV1_min = 0; %0
bound.N_REV1_max = 2; %3
% N REV 2 (7)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 2; %3
% N REV 3 (8)
bound.N_REV3_min = 0; %0
bound.N_REV3_max = 2; %3
% N REV 4 (9)
bound.N_REV4_min = 0; %0
bound.N_REV4_max = 2; %3
% N REV GA (20)
bound.N_REV_GA_min = 0; %0
bound.N_REV_GA_max = 1; %3
% ID Permutation (10)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 
% bound.IDP_max = length(data.asteroid_names);
% C3 stuff
% Constraint on C3 Launcher (11)
sim.C3_max = 30; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (12)
bound.az_min = -pi;
bound.az_max = pi;
% elevation (13)
bound.el_min = -pi/2;
bound.el_max = pi/2;
% GA Stuff in
% v_inf_magn2 (15)
bound.v_inf_magn2_min = 0;
bound.v_inf_magn2_max = 15/sim.DU*sim.TU;
% azimuth2 (16)
bound.az2_min = -pi;
bound.az2_max = pi;
% elevation2 (17)
bound.el2_min = -pi/2;
bound.el2_max = pi/2;
% GA Stuff out
% azimuth3 (18)
bound.az3_min = -pi;
bound.az3_max = pi;
% elevation3 (19)
bound.el3_min = -pi/2;
bound.el3_max = pi/2;

% FB on asteroids
% ASTEROID 1 FLYBY
% vinf_ast1 (21)
bound.vinf_ast1_min = 0;
bound.vinf_ast1_max = 2*sim.TU/sim.DU;
% azimuth_ast1_in (22)
bound.az_ast1_min = -pi;
bound.az_ast1_max = pi;
% elevation_ast1_in (23)
bound.el_ast1_min = -pi/2;
bound.el_ast1_max = pi/2;

% ASTEROID 2 FLYBY
% vinf_ast2 (24)
bound.vinf_ast2_min = 0;
bound.vinf_ast2_max = 2*sim.TU/sim.DU;
% azimuth_ast2_in (25)
bound.az_ast2_min = -pi;
bound.az_ast2_max = pi;
% elevation_ast2_in (26)
bound.el_ast2_min = -pi/2;
bound.el_ast2_max = pi/2;

% ASTEROID 3 FLYBY
% vinf_ast3 (27)
bound.vinf_ast3_min = 0;
bound.vinf_ast3_max = 2*sim.TU/sim.DU;
% azimuth_ast3_in (28)
bound.az_ast3_min = -pi;
bound.az_ast3_max = pi;
% elevation_ast3_in (29)
bound.el_ast3_min = -pi/2;
bound.el_ast3_max = pi/2;

% ASTEROID 4 FLYBY
% vinf_ast4 (30)
bound.vinf_ast4_min = 0;
bound.vinf_ast4_max = 2*sim.TU/sim.DU;
% azimuth_ast1_in (31)
bound.az_ast4_min = -pi;
bound.az_ast4_max = pi;
% elevation_ast1_in (32)
bound.el_ast4_min = -pi/2;
bound.el_ast4_max = pi/2;

% x = [MJD0,TOF1,TOF2,TOF3,TOF4,NREV1,NREV2,NREV3,NREV4,IDP,v_inf_magn,az,el,
%      v_inf_magn2,azimuth2,elevation2,azimuth3,elevation3,NREVGA] and all
%      asteroids flyby variables

bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.TOF3_min, bound.TOF4_min, bound.N_REV1_min, bound.N_REV2_min, ...
    bound.N_REV3_min, bound.N_REV4_min, bound.IDP_min, ...
    bound.v_inf_magn_min, bound.az_min, bound.el_min, bound.TOFGA_min,  ...
    bound.v_inf_magn2_min, bound.az2_min, bound.el2_min, bound.az3_min, ...
    bound.el3_min, bound.N_REV_GA_min, ...
    bound.vinf_ast1_min, bound.az_ast1_min, bound.el_ast1_min, ...
    bound.vinf_ast2_min, bound.az_ast2_min, bound.el_ast2_min, ...
    bound.vinf_ast3_min, bound.az_ast3_min, bound.el_ast3_min, ...
    bound.vinf_ast4_min, bound.az_ast4_min, bound.el_ast4_min]; % Lower bound



bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.TOF3_max, bound.TOF4_max, bound.N_REV1_max, bound.N_REV2_max, ...
    bound.N_REV3_max, bound.N_REV4_max, bound.IDP_max, ...
    bound.v_inf_magn_max, bound.az_max, bound.el_max, bound.TOFGA_max,  ...
    bound.v_inf_magn2_max, bound.az2_max, bound.el2_max, bound.az3_max, ...
    bound.el3_max, bound.N_REV_GA_max, ...
    bound.vinf_ast1_max, bound.az_ast1_max, bound.el_ast1_max, ...
    bound.vinf_ast2_max, bound.az_ast2_max, bound.el_ast2_max, ...
    bound.vinf_ast3_max, bound.az_ast3_max, bound.el_ast3_max, ...
    bound.vinf_ast4_max, bound.az_ast4_max, bound.el_ast4_max]; % Upper bound


%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (6) NREV1, (7) NREV2, (8) NREV3, (9) NREV4, (10) IDP, (20) NREVGA
constr.IntCon = [6,7,8,9,10,20];

%% Options
options = optimoptions(@ga);
% options.PlotFcn = @gaplotpareto;
options.Display = 'iter';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
%options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
% options.HybridFcn = 'fgoalattain';

options.PopulationSize = 1000; % ideal 1000
options.MaxGenerations = 300; % ideal 100

options.FunctionTolerance = 1e-6; %1e-9
options.MaxStallGenerations = ceil(options.MaxGenerations/10);

% Parallel pool
% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;

%% Build the soo
FitnessFunction = @(x) ff_1FL_GA(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

% % Build soo with ps
% options = optimoptions('particleswarm');
% options.HybridFcn = @fmincon;
% options.SwarmSize = 1000; % Default is min(100,10*nvars),
% options.MaxIterations = 300; %  Default is 200*nvars
% options.MaxStallIterations = 50; % Default 20
% options.Display = 'iter';
% options.FunctionTolerance = 1e-6;
% 
% Parallel pool
% Open the parallel pool
% par_pool = gcp; 
% if isempty(par_pool)
%     poolsize = 0;
% else
%     poolsize = par_pool.NumWorkers;
% end
% 
% FitnessFunction = @(x) ff_1FL_GA_Ale_Enri(x,sim,data); % Function handle to the fitness function
% numberOfVariables = length(bound.ub); % Number of decision variables
% [x,Fval,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables,...
%     bound.lb,bound.ub,options)

%% Building the solution structure
sol.asteroid_1 = data.PermutationMatrix(x(10),1);
sol.asteroid_2 = data.PermutationMatrix(x(10),2);
sol.asteroid_3 = data.PermutationMatrix(x(10),3);
sol.asteroid_4 = data.PermutationMatrix(x(10),4);

sol.departure_date = mjd20002date(x(1)*sim.TU/(3600*24));
sol.TOFGA = x(14);
sol.TOF1  = x(2);
sol.TOF2  = x(3);
sol.TOF3  = x(4);
sol.TOF4  = x(5);
sol.mission_duration_d = (sol.TOFGA + sol.TOF1+sol.TOF2+sol.TOF3+sol.TOF4)*sim.TU/(3600*24);
sol.mission_duration_y = sol.mission_duration_d/365;
sol.arrival_date = mjd20002date(x(1)*sim.TU/(3600*24)+sol.mission_duration_d);

sol.objfun = Fval;
sol.Nrev = [x(20), x(6), x(7), x(8), x(9)];
sol.v_inf_magn = x(11)*sim.DU/sim.TU;
sol.az = x(12);
sol.el = x(13);

sol.v_inf_magn2 = x(15)*sim.DU/sim.TU;

sol.v_inf_ast1_magn = x(21)*sim.DU/sim.TU;

sol.v_inf_ast2_magn = x(24)*sim.DU/sim.TU;

sol.v_inf_ast3_magn = x(27)*sim.DU/sim.TU;

sol.v_inf_ast4_magn = x(30)*sim.DU/sim.TU;


%% characteristic quantities plot and Checks
[output, r_encounter, sol] = plot_ff_1FL_GA_Ale(x,sim,data,sol);

figure('Name','Time Evolution of Stuff')
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xline(output.t(sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(2*sim.n_sol)*sim.TU/86400,'LineWidth',2);
% xline(output.t(3*sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(4*sim.n_sol)*sim.TU/86400,'LineWidth',2);
% xline(output.t(5*sim.n_sol)*sim.TU/86400,'LineWidth',2)
% xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.Thrust(:,2));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.Thrust(:,3));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('Mass [kg]')

%% orbit plots
% transfer orbits
r_transf_orbit_GA  = [output.r.legGA.*cos(output.theta.legGA), ...
    output.r.legGA.*sin(output.theta.legGA), output.z.legGA];
R_transf_orbit_GA = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_GA,sim.n_sol,output.Href.legGA);

r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.GA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.ast1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

r_transf_orbit_3  = [output.r.leg3.*cos(output.theta.leg3), ...
    output.r.leg3.*sin(output.theta.leg3), output.z.leg3];
R_transf_orbit_3 = rotate_local2ecplitic(r_encounter.ast2,r_transf_orbit_3,sim.n_sol,output.Href.leg3);

r_transf_orbit_4  = [output.r.leg4.*cos(output.theta.leg4), ...
    output.r.leg4.*sin(output.theta.leg4), output.z.leg4];
R_transf_orbit_4 = rotate_local2ecplitic(r_encounter.ast3,r_transf_orbit_4,sim.n_sol,output.Href.leg4);

figure()
% plot3(R_transf_orbit_GA(:,1),R_transf_orbit_GA(:,2),R_transf_orbit_GA(:,3),...
%     'Color',colors(1,:),'DisplayName','Traj')
 hold on
% hpt1 = plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
%     'Color',colors(1,:));
% hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
%     'Color',colors(1,:));
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% hpt3 = plot3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),...
%     'Color',colors(1,:));
%hpt3.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt4 = plot3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),...
    'Color',colors(1,:));
hpt4.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.GA(1),r_encounter.GA(2),r_encounter.GA(3),...
    '*','Color',colors(2,:),'DisplayName','GA')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.ast2(1),r_encounter.ast2(2),r_encounter.ast2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.ast3(1),r_encounter.ast3(2),r_encounter.ast3(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.ast4(1),r_encounter.ast4(2),r_encounter.ast4(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_4)
% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,2); % mars
% Asteroids
fraction_of_the_orbit = 1/6;
plot_asteorid_orbit(output.t1(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(output.t2(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_2,colors,4);
plot_asteorid_orbit(output.t3(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_3,colors,5);
plot_asteorid_orbit(output.t4(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_4,colors,6);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
legend('show')
view(2)

%% post analysis and checks
max_Thrust = max(output.T_magn)*1000; %[mN]

transfGA_check = abs(sol.TOFGA - output.tGA(end))*sim.TU/(3600*24);
transf1_check = abs(sol.TOF1- output.t1(end))*sim.TU/(3600*24);
transf2_check = abs(sol.TOF2- output.t2(end))*sim.TU/(3600*24);
transfa_check = abs(sol.TOFa - output.ta(end))*sim.TU/(3600*24);
transfb_check = abs(sol.TOFb - output.tb(end))*sim.TU/(3600*24);

Traj = [R_transf_orbit_GA;R_transf_orbit_1;R_transf_orbit_2;R_transf_orbit_3;R_transf_orbit_4];
max_dist_sun = max(vecnorm(Traj,2,2));
min_dist_sun = min(vecnorm(Traj,2,2));

mjd_earth_plot = [output.tGA;output.t1;output.t2;output.t3;output.t4].*sim.TU/86400 + sol.departure_mjd2000;
for k=1:length(mjd_earth_plot)
    [kep,ksun] = uplanet(mjd_earth_plot(k), 3);
    [r__E, ~] = sv_from_coe(kep,ksun);  
    R__E(k,:)=r__E/sim.DU; % it's in km it becomes AU, to be plotted
end
clearvars k r__E

dist_sc_earth = Traj - R__E;
max_dist_earth = max(vecnorm(dist_sc_earth,2,2));
min_dist_earth = min(vecnorm(dist_sc_earth,2,2));

% Propulsion for NIKITA
Propulsion.Magnitude_Thrust = output.T_magn;
Propulsion.InPlane_Thrust   = output.Thrust(:,1);
Propulsion.gamma_angle      = output.Thrust(:,2);
Propulsion.OutofPlane_Thrust = output.Thrust(:,3);