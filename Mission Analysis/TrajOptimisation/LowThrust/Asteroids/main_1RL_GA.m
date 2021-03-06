%% --------------------------------------------------------------------- %%
%% ------------- Earth Ast1 Ast2 Ast3 Ast4 Transfer -------------------- %%
%% ------------------------- ARCH 1RL + GA ----------------------------- %%
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

sim.case_name = 'ARCH 1RL + EARTH/MARS GRAVITY ASSIST';

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

%load('data_processed_9_42_475.mat') 

load('data_processed_42_475.mat')
%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                          % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 200; %100                   % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';    

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse  
sim.M = 100; % SC wet mass [kg]
sim.M_pods = 1; % mass of the payloads+landing stuff [kg]

sim.ID_FLYBY = 3; % flyby planet
sim.c = 1e+3;
sim.max_Available_Thrust = 0.02; % N

AU = astroConstants(2);
%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; %2024
bound.date_ld = [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOFGA (2)
bound.TOFGA_min = 0.4*365*3600*24/sim.TU; %10*3600
bound.TOFGA_max = 2*365*3600*24/sim.TU; 
% TOF1 (3)
bound.TOF1_min = 0.3*365*3600*24/sim.TU; 
bound.TOF1_max = 2*365*3600*24/sim.TU; %2*365
% TOF2 (4)
bound.TOF2_min = 0.3*365*3600*24/sim.TU; 
bound.TOF2_max = 2*365*3600*24/sim.TU; 
% TOF3 (5)
bound.TOF3_min = 0.3*365*3600*24/sim.TU; 
bound.TOF3_max = 2*365*3600*24/sim.TU; 
% TOF4 (6)
bound.TOF4_min = 0.3*365*3600*24/sim.TU; 
bound.TOF4_max = 2*365*3600*24/sim.TU; 
% N REV GA (7)
bound.N_REVGA_min = 0; %0
bound.N_REVGA_max = 1; %3
% N REV 1 (8)
bound.N_REV1_min = 0; %0
bound.N_REV1_max = 1; %3
% N REV 2 (9)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 1; %3
% N REV 3 (10)
bound.N_REV3_min = 0; %0
bound.N_REV3_max = 1; %3
% N REV 4 (11)
bound.N_REV4_min = 0; %0
bound.N_REV4_max = 1; %3
% ID Permutation (12)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 
% Constraint on C3 Launcher (13)
sim.C3_max = 30; % km^2/s^2 %20
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (14)
bound.az_min = -pi;
bound.az_max =  pi;
% elevation (15)
bound.el_min = -pi/2;
bound.el_max =  pi/2;
% vinf2 (GRAVITY ASSIST) (16)
bound.v_inf2_magn_min = 0;
bound.v_inf2_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU; %%%%
% azimuth2 (PRE-GRAVITY ASSIST) (17)
bound.az2_min = -pi;
bound.az2_max =  pi;
% elevation (PRE-GRAVITY ASSIST) (18)
bound.el2_min = -pi/2;
bound.el2_max =  pi/2;
% azimuth (POST-GRAVITY ASSIST)(19)
bound.az3_min = -pi;
bound.az3_max =  pi;
% elevation (POST-GRAVITY ASSIST)(20)
bound.el3_min = -pi/2;
bound.el3_max =  pi/2;
% TOC1 (21) - Coasting time 1
bound.TOC1_min = 30*3600*24/sim.TU; 
bound.TOC1_max = 70*3600*24/sim.TU; %100
% TOC2 (22) - Coasting time 2
bound.TOC2_min = 30*3600*24/sim.TU; 
bound.TOC2_max = 70*3600*24/sim.TU; %100
% TOC3 (23) - Coasting time 3
bound.TOC3_min = 30*3600*24/sim.TU; 
bound.TOC3_max = 70*3600*24/sim.TU; %100


% x = [MJD0,TOFGA, TOF1,TOF2,TOF3,TOF4,NREVGA, NREV1,NREV2,NREV3,NREV4,IDP,v_inf_magn,az,el,v_inf_magn2,az2,el2,az3,el3,TOC1,TOC2,TOC3]
bound.lb = [bound.mjd2000_ed, bound.TOFGA_min, bound.TOF1_min, bound.TOF2_min, ...
    bound.TOF3_min, bound.TOF4_min, bound.N_REVGA_min, bound.N_REV1_min, bound.N_REV2_min, ...
    bound.N_REV3_min, bound.N_REV4_min, bound.IDP_min, ...
    bound.v_inf_magn_min, bound.az_min, bound.el_min, ...
    bound.v_inf2_magn_min, bound.az2_min, bound.el2_min, bound.az3_min, bound.el3_min,bound.TOC1_min,bound.TOC2_min,bound.TOC3_min]; % Lower bound

bound.ub = [bound.mjd2000_ld, bound.TOFGA_max, bound.TOF1_max, bound.TOF2_max, ...
    bound.TOF3_max, bound.TOF4_max, bound.N_REVGA_max, bound.N_REV1_max, bound.N_REV2_max, ...
    bound.N_REV3_max, bound.N_REV4_max, bound.IDP_max, ...
    bound.v_inf_magn_max, bound.az_max, bound.el_max, ...
    bound.v_inf2_magn_max, bound.az2_max, bound.el2_max, bound.az3_max, bound.el3_max,bound.TOC1_max,bound.TOC2_max,bound.TOC3_max ]; % Upper bound


%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% NERVGA (7), NREV1 (8), NREV2 (9), NREV3 (10), NREV4 (11), IDP
constr.IntCon = [7,8,9,10,11,12];

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
FitnessFunction = @(x) ff_1RL_GA(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Building the solution structure
sol.asteroid_1 = data.PermutationMatrix(x(12),1);
sol.asteroid_2 = data.PermutationMatrix(x(12),2);
sol.asteroid_3 = data.PermutationMatrix(x(12),3);
sol.asteroid_4 = data.PermutationMatrix(x(12),4);

sol.departure_mjd2000 = x(1)*sim.TU/(3600*24);
sol.departure_date = mjd20002date(x(1)*sim.TU/(3600*24));
% sol.TOFGA = x(2)*sim.TU/86400;
% sol.TOF1  = x(3)*sim.TU/86400;
% sol.TOF2  = x(4)*sim.TU/86400;
% sol.TOF3  = x(5)*sim.TU/86400;
% sol.TOF4  = x(6)*sim.TU/86400;
sol.TOFGA = x(2);
sol.TOF1  = x(3);
sol.TOF2  = x(4);
sol.TOF3  = x(5);
sol.TOF4  = x(6);
sol.TOC1  = x(21);
sol.TOC2  = x(22);
sol.TOC3  = x(23);

sol.mission_duration_d = (sol.TOFGA + sol.TOF1+ sol.TOC1 + sol.TOF2 + sol.TOC2 + sol.TOF3 + sol.TOC3 + sol.TOF4)*sim.TU/86400;
sol.mission_duration_y = sol.mission_duration_d/365;
sol.arrival_date = mjd20002date(x(1)*sim.TU/(3600*24)+sol.mission_duration_d);
sol.Nrev = [x(7), x(8), x(9), x(10), x(11)];
sol.v_inf_magn = x(13)*sim.DU/sim.TU;
sol.az = x(14);
sol.el = x(15);
sol.v_inf2_magn = x(16)*sim.DU/sim.TU;
sol.az2 = x(17);
sol.el2 = x(18);
sol.az3 = x(19);
sol.el3 = x(20);
sol.obj_fun = Fval;

%% characteristic quantities plot
[output, r_encounter,r_departure,R_coasting, v_encounter, sol] = plot_ff_1RL_GA(x,sim,data, sol);

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
xline(sol.TOFGA*sim.TU/86400,'LineWidth',2);
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
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
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
% xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('Mass [kg]')

%% orbit plots
% transfer orbits
r_transf_orbit_GA = [output.r.GA.*cos(output.theta.GA), ...
    output.r.GA.*sin(output.theta.GA), output.z.GA];
R_transf_orbit_GA = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_GA,sim.n_sol,output.Href.GA);

r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.GA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_departure.ast1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

r_transf_orbit_3  = [output.r.leg3.*cos(output.theta.leg3), ...
    output.r.leg3.*sin(output.theta.leg3), output.z.leg3];
R_transf_orbit_3 = rotate_local2ecplitic(r_departure.ast2,r_transf_orbit_3,sim.n_sol,output.Href.leg3);

r_transf_orbit_4  = [output.r.leg4.*cos(output.theta.leg4), ...
    output.r.leg4.*sin(output.theta.leg4), output.z.leg4];
R_transf_orbit_4 = rotate_local2ecplitic(r_departure.ast3,r_transf_orbit_4,sim.n_sol,output.Href.leg4);

figure()
plot3(R_transf_orbit_GA(:,1),R_transf_orbit_GA(:,2),R_transf_orbit_GA(:,3),...
    'Color',colors(7,:),'DisplayName','Traj GA')
hold on
% plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
%     'Color',colors(1,:),'DisplayName','Traj')
% 
plot3( R_coasting.ast1(:,1), R_coasting.ast1(:,2), R_coasting.ast1(:,3),...
    '--','Color',colors(1,:),'DisplayName','Coast 1')

% hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
%     'Color',colors(1,:));
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3( R_coasting.ast2(:,1), R_coasting.ast2(:,2), R_coasting.ast2(:,3),...
     '--','Color',colors(1,:),'DisplayName','Coast 2')
% % 
plot3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),...
     'Color',colors(1,:));

plot3( R_coasting.ast3(:,1), R_coasting.ast3(:,2), R_coasting.ast3(:,3),...
     '--','Color',colors(1,:),'DisplayName','Coast 3')
% % 
% plot3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),...
%      'Color',colors(1,:));


plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.GA(1),r_encounter.GA(2),r_encounter.GA(3),...
    '*','Color',colors(9,:),'DisplayName','Gravity assist')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_departure.ast1(1),r_departure.ast1(2),r_departure.ast1(3),...
    '*','Color',colors(5,:),'DisplayName','Dep ast1')
plot3(r_encounter.ast2(1),r_encounter.ast2(2),r_encounter.ast2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_departure.ast2(1),r_departure.ast2(2),r_departure.ast2(3),...
    '*','Color',colors(11,:),'DisplayName','Dep ast2')
plot3(r_encounter.ast3(1),r_encounter.ast3(2),r_encounter.ast3(3),...
    '^','Color',colors(12,:),'DisplayName',sol.asteroid_3)
plot3(r_departure.ast3(1),r_departure.ast3(2),r_departure.ast3(3),...
    '*','Color',colors(10,:),'DisplayName','Dep ast3')
plot3(r_encounter.ast4(1),r_encounter.ast4(2),r_encounter.ast4(3),...
    '^','Color',colors(9,:),'DisplayName',sol.asteroid_4)

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 

% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,2); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000 + output.tEnd.Leg12;
hello_orbit2 = sol.departure_mjd2000 + output.tEnd.Leg22;
hello_orbit3 = sol.departure_mjd2000 + output.tEnd.Leg32;
hello_orbit4 = sol.departure_mjd2000 + output.tEnd.Leg42;
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);
plot_asteorid_orbit(hello_orbit3,fraction_of_the_orbit,sol.asteroid_3,colors,5);
plot_asteorid_orbit(hello_orbit4,fraction_of_the_orbit,sol.asteroid_4,colors,6);

% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

% Check position Earth at the end
MJD_end = sol.output_4.t(end)*sim.TU/(3600*24);
[kep_EA_end,ksun] = uplanet(MJD_end, sim.ID_FLYBY);
[rEA_end, vEA_end] = sv_from_coe(kep_EA_end,ksun);

plot3(rEA_end(1)/AU,rEA_end(2)/AU,rEA_end(3)/AU,'r*','DisplayName','Earth end')

%% post analysis and Checks
max_Thrust = max(sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2))*1000; %[mN]

transfGA_check = abs(sol.TOFGA-sol.output_GA.t(end))*sim.TU/(3600*24);
transf1_check = abs(sol.TOF1-sol.output_1.t(end))*sim.TU/(3600*24);
transf2_check = abs(sol.TOF2-sol.output_2.t(end))*sim.TU/(3600*24);
transf3_check = abs(sol.TOF3 -sol.output_3.t(end))*sim.TU/(3600*24);
transf4_check = abs(sol.TOF4 -sol.output_4.t(end))*sim.TU/(3600*24);

Traj = [R_transf_orbit_GA;R_transf_orbit_1;R_transf_orbit_2;R_transf_orbit_3;R_transf_orbit_4];
max_dist_sun = max(vecnorm(Traj,2,2));
min_dist_sun = min(vecnorm(Traj,2,2));

mjd_earth_plot = output.t.*sim.TU/86400 + sol.departure_mjd2000;
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
Propulsion.Magnitude_Thrust = sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2);
Propulsion.InPlane_Thrust   = output.Thrust(:,1);
Propulsion.gamma_angle      = output.Thrust(:,2);
Propulsion.OutofPlane_Thrust = output.Thrust(:,3);

%% Plot with thrust vectors
Tlocal_transf_orbit_GA  = [-sol.output_GA.Thrust(:,1).*sin(output.theta.GA), ...
    sol.output_GA.Thrust(:,1).*cos(output.theta.GA), sol.output_GA.Thrust(:,3)];
T_transf_orbit_GA = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit_GA,sim.n_sol,output.Href.GA);

Tlocal_transf_orbit_1  = [-sol.output_1.Thrust(:,1).*sin(output.theta.leg1), ...
    sol.output_1.Thrust(:,1).*cos(output.theta.leg1), sol.output_1.Thrust(:,3)];
T_transf_orbit_1 = rotate_local2ecplitic(r_encounter.GA,Tlocal_transf_orbit_1,sim.n_sol,output.Href.leg1);


T_coasting_1 = zeros(size(R_coasting.ast1));

Tlocal_transf_orbit_2  = [-sol.output_2.Thrust(:,1).*sin(output.theta.leg2), ...
    sol.output_2.Thrust(:,1).*cos(output.theta.leg2), sol.output_2.Thrust(:,3)];
T_transf_orbit_2 = rotate_local2ecplitic(r_departure.ast1,Tlocal_transf_orbit_2,sim.n_sol,output.Href.leg2);


figure()
plot3(R_transf_orbit_GA(:,1),R_transf_orbit_GA(:,2),R_transf_orbit_GA(:,3),...
    'Color',colors(7,:),'DisplayName','Traj GA')
hold on
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj')
% 
% plot3( R_coasting.ast1(:,1)./AU, R_coasting.ast1(:,2)./AU, R_coasting.ast1(:,3)./AU,...
%     '--','Color',colors(1,:),'DisplayName','Coast')

hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.GA(1),r_encounter.GA(2),r_encounter.GA(3),...
    '*','Color',colors(9,:),'DisplayName','Gravity assist')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_departure.ast1(1),r_departure.ast1(2),r_departure.ast1(3),...
    '^','Color',colors(5,:),'DisplayName','Dep ast1')
plot3(r_encounter.ast2(1),r_encounter.ast2(2),r_encounter.ast2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)

quiver3(R_transf_orbit_GA(:,1),R_transf_orbit_GA(:,2),R_transf_orbit_GA(:,3),T_transf_orbit_GA(:,1),T_transf_orbit_GA(:,2),T_transf_orbit_GA(:,3),2)
hold on
quiver3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3),2)
quiver3(R_coasting.ast1(:,1)./AU,R_coasting.ast1(:,2)./AU,R_coasting.ast1(:,3)./AU,T_coasting_1(:,1),T_coasting_1(:,2),T_coasting_1(:,3),2)
quiver3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3),2)

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
title('In-plane + out-of-plane Thrust')
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
%view(2)

figure;
plot(R_transf_orbit_GA(:,1),R_transf_orbit_GA(:,2),...
    'Color',colors(7,:),'DisplayName','Traj GA')
hold on
plot(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),...
    'Color',colors(1,:),'DisplayName','Traj')

plot( R_coasting.ast1(:,1)./AU, R_coasting.ast1(:,2)./AU,...
    '--','Color',colors(1,:),'DisplayName','Coast')

hpt2 = plot(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(r_encounter.EA(1),r_encounter.EA(2),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot(r_encounter.GA(1),r_encounter.GA(2),...
    '*','Color',colors(9,:),'DisplayName','Gravity assist')
plot(r_encounter.ast1(1),r_encounter.ast1(2),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot(r_departure.ast1(1),r_departure.ast1(2),...
    '^','Color',colors(5,:),'DisplayName','Dep ast1')
plot(r_encounter.ast2(1),r_encounter.ast2(2),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
quiver(R_transf_orbit_GA(:,1),R_transf_orbit_GA(:,2),T_transf_orbit_GA(:,1),T_transf_orbit_GA(:,2),2)
hold on
quiver(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),2)
quiver(R_coasting.ast1(:,1)./AU,R_coasting.ast1(:,2)./AU,T_coasting_1(:,1),T_coasting_1(:,2),2)
quiver(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),2)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]');
title('In-plane Thrust(shall be tangential)')
% Sun
plot(0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

Thrust_Helio = [T_transf_orbit_GA(:,1),T_transf_orbit_GA(:,2),T_transf_orbit_GA(:,3);...
    T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3);...
    T_coasting_1(:,1),T_coasting_1(:,2),T_coasting_1(:,3);...
    T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3)];
    