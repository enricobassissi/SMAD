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
load('data_elements_matrix.mat')
[data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
            sequences_local_pruning(data_elements_matrix);

[data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

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
% Departure dates
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1
bound.TOF1_min = 550*3600*24/sim.TU; %600
bound.TOF1_max = 900*3600*24/sim.TU; 
% TOF2
bound.TOF2_min = 1*365*3600*24/sim.TU; %600
bound.TOF2_max = 3*365*3600*24/sim.TU; 
% N REV
bound.N_REV_min = 1; %0
bound.N_REV_max = 1; %3
% ID Permutation
bound.IDP_min = 1; 
% bound.IDP_max = data.HowMany; 
bound.IDP_max = length(data.asteroid_names);
% C3 stuff
% Constraint on C3 Launcher
sim.C3_max = 10; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
bound.alpha_min = -pi;
bound.alpha_max = pi;
bound.beta_min = -pi;
bound.beta_max = pi;

% x = [MJD0,TOF1,TOF2,NREV,IDP,v_inf_magn,alpha,beta]
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.N_REV_min, bound.IDP_min, bound.v_inf_magn_min, ...
    bound.alpha_min, bound.beta_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.N_REV_max, bound.IDP_max, bound.v_inf_magn_max, ...
    bound.alpha_max, bound.beta_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
constr.IntCon = [4,5];

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
FitnessFunction = @(x) ff_ea_ast_LT_soo_NLI(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

dep_opt = mjd20002date(x(1)*sim.TU/(3600*24))
asteroid_1 = data.asteroid_names(x(5))

%% plot
[output, r1_true, r2_true,v1_true,v2_true] = plot_ff_ea_ast_LT_soo_NLI(x,sim,data);

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.Thrust(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.Thrust(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

%%
%JD_departure = x(knee_sol,1);
day1 = [2028 1 1 0 0 0];
day2 = [2031 1 1 0 0 0];

t1 = date2mjd2000(day1);
t2 = date2mjd2000(day2);
times = linspace(t1,t2,1000);

for i=1:length(times)
    % Orbit 1
    [kep1,ksun] = uplanet(times(i),3);
    [r1(i,1:3),v1(i,1:3)] = sv_from_coe(kep1,ksun);
    r1(i,1:3) = r1(i,1:3)/sim.DU;
    
    % Orbit 2
    [kep2,ksun] = uplanet(times(i),4);
    [r2(i,1:3),~] = sv_from_coe(kep2,ksun);
    r2(i,1:3) = r2(i,1:3)/sim.DU;
    
    % Orbit 3
    [kep_ast_1] = uNEO2(times(i),asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [r_ast1(i,1:3), v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    r_ast1(i,1:3) = r_ast1(i,1:3)/sim.DU;
end

r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
R3 = rotate_local2ecplitic(r1_true,r3,sim.n_sol,output.Href);

figure()
plot3(R3(:,1),R3(:,2),R3(:,3),'DisplayName','Traj')
hold on
plot3(r1_true(1),r1_true(2),r1_true(3),'*m','DisplayName','Dep')
plot3(r2_true(1),r2_true(2),r2_true(3),'*c','DisplayName','Arr')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% Earth
plot3(r1(:,1),r1(:,2),r1(:,3),'--g','DisplayName','Earth'); % geocentric equatorial frame 
% Mars
plot3(r2(:,1),r2(:,2),r2(:,3),'--r','DisplayName','Mars');
% Asteroid 1
plot3(r_ast1(:,1),r_ast1(:,2),r_ast1(:,3),'--b','DisplayName','Ast1');
% Sun
plot3(0,0,0,'oy','DisplayName','Sun')
legend('show')

