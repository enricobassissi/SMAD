%% --------------------------------------------------------------------- %%
%% ------------------------- Ast1 Ast2 Transfer ------------------------ %%
%% ------------------------- ARCH 1+4, LT FLYBY ------------------------ %%
%% ------------------------------ SOO ---------------------------------- %%
%% --------------- ONE ASTEROID TO ANOTHER, BEST CASE ------------------ %%
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
% load('data_elements_matrix_9.mat')
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix);
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

load('data_processed_9_1.mat')

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

sim.PS.Isp = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 100; % SC mass [kg]

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 0.5*365*3600*24/sim.TU; 
bound.TOF1_max = 3*365*3600*24/sim.TU; 
% N REV (3)
bound.N_REV_min = 0; %0
bound.N_REV_max = 1; %3
% % C3 stuff
% % Constraint on C3 Launcher (4)
% sim.C3_max = 10; % km^2/s^2
% bound.v_inf_magn_min = 0;
% bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% % azimuth (5)
% bound.az_min = -pi;
% bound.az_max = pi;
% % elevation (6)
% bound.el_min = -pi/2;
% bound.el_max = pi/2;
% ID Permutation (4)
bound.IDP_min = 1; 
% bound.IDP_max = data.HowMany; 
bound.IDP_max = length(data.asteroid_names);
% ID Permutation2 (5)
bound.IDP2_min = 1; 
% bound.IDP_max = data.HowMany; 
bound.IDP2_max = length(data.asteroid_names)-1;

% x = [MJD0,TOF1,NREV,v_inf_magn,az,el,IDP,IDP2]
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, ...
    bound.N_REV_min, bound.IDP_min,bound.IDP2_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, ...
    bound.N_REV_max, bound.IDP_max,bound.IDP2_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
constr.IntCon = [3,4,5];

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

options.PopulationSize = 600; % ideal 1000
options.MaxGenerations = 30; % ideal 100

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
FitnessFunction = @(x) ff_1ast_to_another_ast_LT_soo_NLI(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

dep_opt = mjd20002date(x(1)*sim.TU/(3600*24))
sol.asteroid_1 = data.asteroid_names(x(4))
sol.asteroid_2 = data.asteroid_names(x(5))

%% plot
[output, r_encounter, v_encounter] = plot_ff_1ast_to_another_ast_LT_soo_NLI(x,sim,data);

figure()
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

%% plot orbits
% transfer orbits
r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
R3 = rotate_local2ecplitic(r_encounter.ast1,r3,sim.n_sol,output.Href);

figure()
plot3(R3(:,1),R3(:,2),R3(:,3),'DisplayName','Traj')
hold on
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.ast2(1),r_encounter.ast2(2),r_encounter.ast2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 

% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars
% Asteroids
fraction_of_the_orbit = 1;
plot_asteorid_orbit(output.t(1)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(output.t(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_2,colors,4);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)