%% --------------------------------------------------------------------- %%
%% ------------------------- Earth Ast1 Transfer ----------------------- %%
%% ------------------------- ARCH 1+4, LT FLYBY ------------------------ %%
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

sim.case_name = 'ARCH ID 6: LOW THRUST FLYBY ON EVERY ASTEROID';

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
% load('data_elements_matrix_42.mat')
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix);
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

load('data_processed_42_475.mat')
load('power_propulsion_data.mat')

%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 200;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                     % direction of integration (1 FW, -1 BW)

% sim.vinf = 0; % Parabolic escape

sim.TOF_imposed_flag = 1;

sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
sim.M = 100; % SC mass [kg]
% sim.M = 80; % SC mass [kg]
sim.max_Available_Thrust = 0.05; % [N]

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 0.2*365*3600*24/sim.TU; %600
bound.TOF1_max = 3*365*3600*24/sim.TU; 
% N REV (3)
bound.N_REV_min = 0; %0
bound.N_REV_max = 2; %3
% C3 stuff
% Constraint on C3 Launcher (4)
sim.C3_max = 30; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (5)
bound.az_min = -pi;
bound.az_max = pi;
% elevation (6)
bound.el_min = -pi/2;
bound.el_max = pi/2;
% ID Permutation (7)
bound.IDP_min = 1; 
% bound.IDP_max = data.HowMany; 
bound.IDP_max = length(data.asteroid_names);

% x = [MJD0,TOF1,NREV,v_inf_magn,az,el,IDP]
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, ...
    bound.N_REV_min, bound.v_inf_magn_min, ...
    bound.az_min, bound.el_min, bound.IDP_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, ...
    bound.N_REV_max, bound.v_inf_magn_max, ...
    bound.az_max, bound.el_max, bound.IDP_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
constr.IntCon = [3,7];

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
options.MaxGenerations = 50; % ideal 100

options.FunctionTolerance = 1e-6; %1e-9
options.MaxStallGenerations = ceil(options.MaxGenerations/5);

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
FitnessFunction = @(x) ff_ea_1ast_LT_soo_NLI_melia(x,sim,data,power_propulsion_data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Solution Struct
sol.departure_mjd2000 = x(1)*sim.TU/(3600*24);
sol.dep_opt = mjd20002date(x(1)*sim.TU/(3600*24));
sol.TOF = x(2)*sim.TU/(3600*24);
sol.asteroid_1 = data.asteroid_names(x(7));

%% plot
[output,r_encounter,v_encounter,sol] = plot_ff_ea_1ast_LT_soo_NLI_melia(x,sim,data,power_propulsion_data,sol);

figure('Name','Params')
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
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
xlabel('Time [days]')
ylabel('Mass [kg]')

%% plot thrust melia
figure('Name','Thrust Available and Requested')
plot(sol.t,sol.T_available,'DisplayName','T Available')
hold on
plot(sol.t,sol.T_magn_Helio,'DisplayName','T Required')
legend('show')

%% plot orbit
r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
R3 = rotate_local2ecplitic(r_encounter.EA,r3,sim.n_sol,output.Href);

figure('Name','Orbits')
plot3(R3(:,1),R3(:,2),R3(:,3),'DisplayName','Traj')
hold on
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),'*m','DisplayName','Dep')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),'*c','DisplayName','Arr')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(sol.departure_mjd2000,3,colors,8); % earth
plot_planet_orbit(sol.departure_mjd2000,4,colors,2); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000 + output.t(end);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

%% post analysis and checks
max_Thrust = max(sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2))*1000; %[mN]

transf_check = abs(sol.TOF - output.t(end)*sim.TU/(3600*24));

Traj = R3;
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
Tlocal_transf_orbit  = [-output.Thrust(:,1).*sin(output.theta), ...
    output.Thrust(:,1).*cos(output.theta), output.Thrust(:,3)];
Thrust_Helio = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit,sim.n_sol,output.Href);

figure('Name','Thrust Vector Plot')
plot3(R3(:,1),R3(:,2),R3(:,3),...
    'Color',colors(1,:),'DisplayName','Traj')
hold on

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)

quiver3(R3(:,1),R3(:,2),R3(:,3),Thrust_Helio(:,1),Thrust_Helio(:,2),Thrust_Helio(:,3),2)

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
title('In-plane + out-of-plane Thrust')
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
%view(2)

figure('Name','Plane thrust')
hold on
plot(R3(:,1),R3(:,2),...
    'Color',colors(1,:),'DisplayName','Traj')

plot(r_encounter.EA(1),r_encounter.EA(2),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot(r_encounter.ast1(1),r_encounter.ast1(2),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)

quiver(R3(:,1),R3(:,2),Thrust_Helio(:,1),Thrust_Helio(:,2),2,...
    'Color',colors(3,:),'DisplayName',sol.asteroid_1)

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]');
title('In-plane Thrust(shall be tangential)')
% Sun
plot(0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

