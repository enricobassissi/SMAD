%% --------------------------------------------------------------------- %%
%% ---------------------- Earth Ast1 Ast2 Transfer --------------------- %%
%% ------------------------- ARCH 1+4, LT RV --------------------------- %%
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

sim.case_name = 'ARCH ID 6: LOW THRUST RENDEZVOUS ON EVERY ASTEROID + COASTING TIME';

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

% load('data_processed_9_4.mat')
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
sim.PS.Isp = 3450/sim.TU;  % non-dimensional specific impulse
sim.M = 300; % SC wet mass [kg]
sim.M_pods = 5; % mass of the payloads+landing stuff [kg]
sim.max_Available_Thrust = 0.06; % 10 mN, BepiColombo was 250 mN but was much bigger

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 0.3*3600*24/sim.TU; 
bound.TOF1_max = 3*365*3600*24/sim.TU; 
% TOF2 (3)
bound.TOF2_min = 0.3*365*3600*24/sim.TU; 
bound.TOF2_max = 3*365*3600*24/sim.TU; 
% N REV (4)
bound.N_REV1_min = 0; %0
bound.N_REV1_max = 2; %3
% N REV 2 (5)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 2; %3
% ID Permutation (6)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 
% bound.IDP_max = length(data.asteroid_names);
% C3 stuff
% Constraint on C3 Launcher (7)
sim.C3_max = 30; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (8)
bound.alpha_min = -deg2rad(180);
bound.alpha_max = deg2rad(180);
% elvation (9)
bound.beta_min = -deg2rad(90);
bound.beta_max = deg2rad(90);
% coasting time1 (10)
bound.CT1_min = 30*3600*24/sim.TU;
bound.CT1_max = 80*3600*24/sim.TU;
% coasting time2 (11)
bound.CT2_min = 30*3600*24/sim.TU;
bound.CT2_max = 80*3600*24/sim.TU;
% TOF3 (12)
bound.TOF3_min = 0.3*365*3600*24/sim.TU; %600
bound.TOF3_max = 3*365*3600*24/sim.TU; 
% N REV 3 (13)
bound.N_REV3_min = 0; %0
bound.N_REV3_max = 2; %3
% coasting time3 (14)
bound.CT3_min = 30*3600*24/sim.TU;
bound.CT3_max = 80*3600*24/sim.TU;
% TOF4 (15)
bound.TOF4_min = 0.3*365*3600*24/sim.TU; %600
bound.TOF4_max = 3*365*3600*24/sim.TU; 
% N REV 4 (16)
bound.N_REV4_min = 0; %0
bound.N_REV4_max = 2; %3

% x = [MJD0,TOF1,TOF2,NREV,NREV2,IDP,v_inf_magn,alpha,beta,CT1,CT2,TOF3,NREV3,
%      CT3,TOF4,NREV4]
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.N_REV1_min, bound.N_REV2_min, bound.IDP_min, ...
    bound.v_inf_magn_min, bound.alpha_min, bound.beta_min, ...
    bound.CT1_min, bound.CT2_min, bound.TOF3_min, bound.N_REV3_min, ...
    bound.CT3_min, bound.TOF4_min, bound.N_REV4_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.N_REV1_max, bound.N_REV2_max, bound.IDP_max, ...
    bound.v_inf_magn_max, bound.alpha_max, bound.beta_max, ...
    bound.CT1_max, bound.CT2_max, bound.TOF3_max, bound.N_REV3_max, ...
    bound.CT3_max, bound.TOF4_max, bound.N_REV4_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (4) NREV1, (5) NREV2, (6) IDP, (13) NREV3, (16) NREV4 
constr.IntCon = [4,5,6,13,16];

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

options.PopulationSize = 1200; % ideal 1000
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
FitnessFunction = @(x) ff_1RL(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Building the solution structure
sol.asteroid_1 = data.PermutationMatrix(x(6),1);
sol.asteroid_2 = data.PermutationMatrix(x(6),2);
sol.asteroid_3 = data.PermutationMatrix(x(6),3);
sol.asteroid_4 = data.PermutationMatrix(x(6),4);

sol.departure_mjd2000 = x(1)*sim.TU/(3600*24);
sol.departure_date = mjd20002date(x(1)*sim.TU/(3600*24));
sol.TOF1 = x(2)*sim.TU/86400;
sol.TOF2 = x(3)*sim.TU/86400;
sol.TOF3 = x(12)*sim.TU/86400;
sol.TOF4 = x(15)*sim.TU/86400;
sol.CT1 = x(10)*sim.TU/86400;
sol.CT2 = x(11)*sim.TU/86400;
sol.CT3 = x(14)*sim.TU/86400;
sol.mission_duration_d = sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3+sol.CT3+sol.TOF4;
sol.mission_duration_y = sol.mission_duration_d/365;
sol.arrival_date = mjd20002date(x(1)*sim.TU/(3600*24)+sol.mission_duration_d);

sol.obj_fun = Fval;
sol.Nrev = [x(4), x(5), x(13), x(16)];
sol.v_inf_magn = x(7)*sim.DU/sim.TU;
sol.az = x(8);
sol.el = x(9);

%% characteristic quantities plot
[output, r_encounter, v_encounter, sol] = plot_ff_1RL(x,sim,data,sol);

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
xline(output.tEnd.Leg11,'LineWidth',2); xline(output.tEnd.Leg12,'LineWidth',2);
xline(output.tEnd.Leg21,'LineWidth',2); xline(output.tEnd.Leg22,'LineWidth',2);
xline(output.tEnd.Leg31,'LineWidth',2); xline(output.tEnd.Leg32,'LineWidth',2);
xline(output.tEnd.Leg41,'LineWidth',2); xline(output.tEnd.Leg42,'LineWidth',2);
% xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.Thrust(:,2));
% xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.Thrust(:,3));
% xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
% xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
xline(output.tEnd.Leg11,'LineWidth',2); xline(output.tEnd.Leg12,'LineWidth',2);
xline(output.tEnd.Leg21,'LineWidth',2); xline(output.tEnd.Leg22,'LineWidth',2);
xline(output.tEnd.Leg31,'LineWidth',2); xline(output.tEnd.Leg32,'LineWidth',2);
xline(output.tEnd.Leg41,'LineWidth',2); xline(output.tEnd.Leg42,'LineWidth',2);
% xlabel('Time [days]')
ylabel('Mass [kg]')

%% orbit plots
% transfer orbits
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

r_transf_orbit_3  = [output.r.leg3.*cos(output.theta.leg3), ...
    output.r.leg3.*sin(output.theta.leg3), output.z.leg3];
R_transf_orbit_3 = rotate_local2ecplitic(r_encounter.astD2,r_transf_orbit_3,sim.n_sol,output.Href.leg3);

r_transf_orbit_4  = [output.r.leg4.*cos(output.theta.leg4), ...
    output.r.leg4.*sin(output.theta.leg4), output.z.leg4];
R_transf_orbit_4 = rotate_local2ecplitic(r_encounter.astD3,r_transf_orbit_4,sim.n_sol,output.Href.leg4);

figure()
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj')
hold on
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt3 = plot3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),...
    'Color',colors(1,:));
hpt3.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt4 = plot3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),...
    'Color',colors(1,:));
hpt4.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astD2(1),r_encounter.astD2(2),r_encounter.astD2(3),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astA3(1),r_encounter.astA3(2),r_encounter.astA3(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astD3(1),r_encounter.astD3(2),r_encounter.astD3(3),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astA4(1),r_encounter.astA4(2),r_encounter.astA4(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_4)
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

%% post analysis and checks
max_Thrust = max(sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2))*1000; %[mN]

transf1_check = abs(sol.TOF1 - output.t1(end)*sim.TU/(3600*24));
transf2_check = abs(sol.TOF2 - output.t2(end)*sim.TU/(3600*24));
transf3_check = abs(sol.TOF3 - output.t3(end)*sim.TU/(3600*24));
transf4_check = abs(sol.TOF4 - output.t4(end)*sim.TU/(3600*24));

Traj = [R_transf_orbit_1;R_transf_orbit_2;R_transf_orbit_3;R_transf_orbit_4];
max_dist_sun = max(vecnorm(Traj,2,2));
min_dist_sun = min(vecnorm(Traj,2,2));

mjd_earth_plot = [output.t1;output.t2;output.t3;output.t4].*sim.TU/86400 + sol.departure_mjd2000;
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
Tlocal_transf_orbit_1  = [-sol.T_1(:,1).*sin(output.theta.leg1), ...
    sol.T_1(:,1).*cos(output.theta.leg1), sol.T_1(:,3)];
T_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit_1,sim.n_sol,output.Href.leg1);

Tlocal_transf_orbit_2  = [-sol.T_2(:,1).*sin(output.theta.leg2), ...
    sol.T_2(:,1).*cos(output.theta.leg2), sol.T_2(:,3)];
T_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,Tlocal_transf_orbit_2,sim.n_sol,output.Href.leg2);

Tlocal_transf_orbit_3  = [-sol.T_3(:,1).*sin(output.theta.leg3), ...
    sol.T_3(:,1).*cos(output.theta.leg3), sol.T_3(:,3)];
T_transf_orbit_3 = rotate_local2ecplitic(r_encounter.astD2,Tlocal_transf_orbit_3,sim.n_sol,output.Href.leg3);

Tlocal_transf_orbit_4  = [-sol.T_4(:,1).*sin(output.theta.leg4), ...
    sol.T_4(:,1).*cos(output.theta.leg4), sol.T_4(:,3)];
T_transf_orbit_4 = rotate_local2ecplitic(r_encounter.astD3,Tlocal_transf_orbit_4,sim.n_sol,output.Href.leg4);

figure('Name','Thrust Plot')
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj')
hold on
hpt1 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt2 = plot3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt3 = plot3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),...
    'Color',colors(1,:));
hpt3.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astD2(1),r_encounter.astD2(2),r_encounter.astD2(3),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astA3(1),r_encounter.astA3(2),r_encounter.astA3(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astD3(1),r_encounter.astD3(2),r_encounter.astD3(3),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astA4(1),r_encounter.astA4(2),r_encounter.astA4(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_4)

quiver3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3),2)
quiver3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3),2)
quiver3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),T_transf_orbit_3(:,1),T_transf_orbit_3(:,2),T_transf_orbit_3(:,3),2)
quiver3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),T_transf_orbit_4(:,1),T_transf_orbit_4(:,2),T_transf_orbit_4(:,3),2)

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
title('In-plane + out-of-plane Thrust')
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
%view(2)

figure('Name','Plane thrust')
hold on
% plot(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),...
%     'Color',colors(1,:),'DisplayName','Traj')
% hpt1 = plot(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),...
%     'Color',colors(1,:));
% hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),...
%     'Color',colors(1,:),'DisplayName','Traj')
hpt2 = plot(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),...
    'Color',colors(1,:));
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(r_encounter.EA(1),r_encounter.EA(2),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot(r_encounter.astA1(1),r_encounter.astA1(2),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot(r_encounter.astD1(1),r_encounter.astD1(2),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot(r_encounter.astA2(1),r_encounter.astA2(2),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_2)
plot(r_encounter.astD2(1),r_encounter.astD2(2),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_2)
plot(r_encounter.astA3(1),r_encounter.astA3(2),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_3)
plot(r_encounter.astD3(1),r_encounter.astD3(2),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_3)
plot(r_encounter.astA4(1),r_encounter.astA4(2),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_4)

% quiver(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),2,...
%     'Color',colors(3,:),'DisplayName',sol.asteroid_1)
% quiver(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),2,...
%     'Color',colors(4,:),'DisplayName',sol.asteroid_2)
% quiver(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),T_transf_orbit_3(:,1),T_transf_orbit_3(:,2),2,...
%     'Color',colors(5,:),'DisplayName',sol.asteroid_3)
quiver(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),T_transf_orbit_4(:,1),T_transf_orbit_4(:,2),2,...
    'Color',colors(6,:),'DisplayName',sol.asteroid_4)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]');
title('In-plane Thrust(shall be tangential)')
% Sun
plot(0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

Thrust_Helio = [T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3);...
                    T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3);...
                    T_transf_orbit_3(:,1),T_transf_orbit_3(:,2),T_transf_orbit_3(:,3);...
                    T_transf_orbit_4(:,1),T_transf_orbit_4(:,2),T_transf_orbit_4(:,3)];