%% --------------------------------------------------------------------- %%
%% -------------- Ma se io faccio un'ottimizzazione a giro ------------- %%
%% --------- escludendo dalla run quello che ho appena visitato -------- %%
%% --------------------------- che succede? ---------------------------- %%
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
sim.max_Available_Thrust = 0.01; % [N]
sim.massPods = 5; % kg

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
sim.asteroid_to_fish = data.asteroid_names;
FitnessFunction = @(x) ff_chesuccseottimizzunoallavolta(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x1,Fval1,exitFlag1,Output1] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp1 = toc/60;

[output1, r_encounter, v_encounter,sol] = plot_ff_chesuccseottimizzunoallavolta(x1,sim,data,sol);

%% and then let's go to the second
sol.time_arr_ast_ADIM = x1(1)+x1(2);
sol.asteroid_1 = data.asteroid_names(x1(7));
sim.asteroid_to_fish = data.asteroid_names(~contains(data.asteroid_names,sol.asteroid_1));
sim.M = output1.m(end) - sim.massPods; % SC mass [kg]
FitnessFunction = @(x) ff_chesucc1_1ast_to_another_ast(x,sim,data,sol);

%% Boundaries 2
% CT1 (1)
bound2.CT1_min = 30*3600*24/sim.TU; %600
bound2.CT1_max = 80*3600*24/sim.TU; 
% TOF1 (2)
bound2.TOF1_min = 0.2*365*3600*24/sim.TU; %600
bound2.TOF1_max = 3*365*3600*24/sim.TU; 
% N REV (3)
bound2.N_REV_min = 0; %0
bound2.N_REV_max = 2; %3
% ID Permutation (4)
bound2.IDP_min = 1; 
% bound.IDP_max = data.HowMany; 
bound2.IDP_max = length(sim.asteroid_to_fish);

% x = [MJD0,TOF1,NREV,v_inf_magn,az,el,IDP]
bound2.lb = [bound2.CT1_min, bound2.TOF1_min, ...
    bound2.N_REV_min, bound2.IDP_min]; % Lower bound
bound2.ub = [bound2.CT1_max, bound2.TOF1_max, ...
    bound2.N_REV_max, bound2.IDP_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
constr.IntCon = [3,4];

%%
numberOfVariables2 = length(bound2.ub);
tic
[x2,Fval2,exitFlag2,Output2] = ga(FitnessFunction,numberOfVariables2,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound2.lb,...
    bound2.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp2 = toc/60;

[output2, r_encounter2, v_encounter2, sol] = plot_ff_chesucc1_1ast_to_another_ast(x2,sim,data,sol);

%% plot
r_1  = [output1.r.*cos(output1.theta), ...
    output1.r.*sin(output1.theta), output1.z];
R_1 = rotate_local2ecplitic(r_encounter.EA,r_1,sim.n_sol,output1.Href);

r_2  = [output2.r.*cos(output2.theta), ...
    output2.r.*sin(output2.theta), output2.z];
R_2 = rotate_local2ecplitic(r_encounter2.ast1,r_2,sim.n_sol,output2.Href);

plot3(R_1(:,1),R_1(:,2),R_1(:,3));
hold on
plot3(R_2(:,1),R_2(:,2),R_2(:,3));
axis equal
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),'*',...
    'DisplayName','Dep Earth')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),'^',...
    'DisplayName','Arr Ast1')
plot3(r_encounter2.ast1(1),r_encounter2.ast1(2),r_encounter2.ast1(3),'*',...
    'DisplayName','Dep Ast1')
plot3(r_encounter2.ast2(1),r_encounter2.ast2(2),r_encounter2.ast2(3),'^',...
    'DisplayName','Arr Ast2')
legend('show')