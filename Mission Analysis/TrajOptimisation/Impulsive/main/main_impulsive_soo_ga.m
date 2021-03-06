%% --------------------------------------------------------------------- %%
%% ------------------ EA Ast1 Ast2 Ast3 Ast4 Transfer ------------------ %%
%% --------------- Impulsive SOO Genetic Algorithm --------------------- %%
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

sim.case_name = 'ARCH ID 1: IMPULSIVE DOUBLE RENDEZVOUS ON EACH ASTEROID';

%% add path of functions and python stuff
path_str=split(pwd, 'TrajOptimisation\Impulsive\main');
path_utils=string(path_str(1))+'Utils';
addpath(genpath(path_utils));
path_py=string(path_str(1))+'PyInterface\NEO_API_py';
addpath(genpath(path_py));
path_neoeph=string(path_str(1))+'NeoEph';
addpath(genpath(path_neoeph));
path_str=split(pwd, 'main');
path_imp=string(path_str(1));
addpath(genpath(path_imp));

%% INTRO ADIMENSIONALISATION
% sim.mu = 1.32712440017987e11; % Sun planetary constant (mu = mass * G) (from DE405) [km^3/s^2]
% sim.DU = 149597870.691; % Distance Unit = Astronomical Unit (AU) (from DE405) [km]
% sim.TU = (sim.DU^3/sim.mu)^0.5; % Time Unit
% sim.mu = 1;

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(path_py+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%% Asteroids
AU = astroConstants(2);
muSun = astroConstants(4);

% data extraction section
data.asteroid_names = ["2006HX57";"2008XU2";"2008KN11";"2012SY49";"2012QD8";"2020UE";...
                  "2006SC";"2005WG57";"2012BY1"];

% Number of possible combination of 4 asteroids among the ones in the list
data.HowMany = factorial(length(data.asteroid_names)) / factorial(length(data.asteroid_names) - 4);
[data.PermutationMatrix, ~] = permnUnique(data.asteroid_names, 4);

%% uNEO
% try 
    load('data.mat')
% catch
    % if the asteroid have changed, run the find_eph_neo below, it takes about 1 min
%     [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);
%     save('data.mat', data);
% end

%% Boundaries
% Departure dates (1)
sim.soo_lim.date_ed = [2022, 1, 1, 0, 0, 0];
sim.soo_lim.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.soo_lim.mjd2000_ed = date2mjd2000(sim.soo_lim.date_ed);
sim.soo_lim.mjd2000_ld = date2mjd2000(sim.soo_lim.date_ld);
% TOF1 (2)
sim.soo_lim.TOF1_min = 200; % days
sim.soo_lim.TOF1_max = 3*365; % days
% Buffer time 1 (3)
sim.soo_lim.bt1_min = 30;
sim.soo_lim.bt1_max = 180;
% TOF2 (4)
sim.soo_lim.TOF2_min = 50; % days
sim.soo_lim.TOF2_max = 3*365; % days
% Matrix of permutations (5)
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.soo_lim.permutations_low = -0.5; 
sim.soo_lim.permutations_up = data.HowMany + 0.4999;
% Buffer time 2 (6)
sim.soo_lim.bt2_min = 30;
sim.soo_lim.bt2_max = 180;
% TOF3 (7)
sim.soo_lim.TOF3_min = 50; % days
sim.soo_lim.TOF3_max = 3*365; % days
% Buffer time 3 (8)
sim.soo_lim.bt3_min = 30;
sim.soo_lim.bt3_max = 180;
% TOF4 (9)
sim.soo_lim.TOF4_min = 50; % days
sim.soo_lim.TOF4_max = 3*365; % days

% x = [MJD0,TOF1,buffer_time,TOF2,ID_permutation,...
%      buffer_time2,TOF3,buffer_time3,TOF4]
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed, sim.soo_lim.TOF1_min,...
      sim.soo_lim.bt1_min,...
      sim.soo_lim.TOF2_min,sim.soo_lim.permutations_low,sim.soo_lim.bt2_min,...
      sim.soo_lim.TOF3_min,sim.soo_lim.bt3_min,sim.soo_lim.TOF4_min]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max,...
      sim.soo_lim.bt1_max,...
      sim.soo_lim.TOF2_max,sim.soo_lim.permutations_up,sim.soo_lim.bt2_max,...
      sim.soo_lim.TOF3_max,sim.soo_lim.bt3_max,sim.soo_lim.TOF4_max]; % Upper bound

% Constraint on C3 Launcher
sim.C3_max = 30; % km^2/s^2

%% Constraints
sim.soo_constr.A = []; % linear inequality constraints
sim.soo_constr.b = []; % linear inequality constraints
sim.soo_constr.Aeq = []; % linear equality constraints
sim.soo_constr.beq = []; % linear equality constraints
sim.soo_constr.nonlcon = []; % linear equality constraints
% sim.soo_constr.IntCon = [5];

%% Options
options = optimoptions(@ga);
% options.PlotFcn = {@gaplotbestf};
options.Display = 'iter';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
% options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
% options.HybridFcn = @fmincon;

options.PopulationSize = 1000; 
options.MaxGenerations = 200; 
options.FunctionTolerance = 1e-6;
options.MaxStallGenerations = 30;

% Parallel pool
% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;
% options.UseParallel = false;
%% Build the soo
FitnessFunction = @(x) ff_impulsive_soo(x, data, sim); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,sim.soo_constr.A, ...
    sim.soo_constr.b,sim.soo_constr.Aeq,sim.soo_constr.beq,sim.soo_bound.lb,sim.soo_bound.ub,...
    sim.soo_constr.nonlcon,options);
el_time_min_pp = toc/60;

%% Build solution structure
% set the knee as main solution
asteroid_sequence = data.PermutationMatrix(round(x(5)),:);
sol.ast_1 = asteroid_sequence(1);
sol.ast_2 = asteroid_sequence(2);
sol.ast_3 = asteroid_sequence(3);
sol.ast_4 = asteroid_sequence(4);
sol.MJD0 = x(1);
sol.dep_date = mjd20002date(sol.MJD0)';
sol.TOF_tot_D = x(2)+x(3)+x(4)+x(6)+x(7)+x(8)+x(9);
sol.TOF_tot_Y = sol.TOF_tot_D/365;
sol.end_of_mission_date = mjd20002date(sol.MJD0+sol.TOF_tot_D)';
sol.dV_tot = Fval(1);
sol.TOF1 = x(2);
sol.buffer_time1 = x(3);
sol.TOF2 = x(4);
sol.buffer_time2 = x(6);
sol.TOF3 = x(7);
sol.buffer_time3 = x(8);
sol.TOF4 = x(9);

%% Mass Consumption for High Thrust Impulsive Case
g0 = 9.81; %m/s^2
% https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
Isp = 230; %s 
m_dry = 100; %kg
m_prop = m_dry*(exp(sol.dV_tot*1e3/(g0*Isp)) - 1); %kg
%% Plot trajectories
sol = plot_mission_4neo_rendezvous(sol,asteroid_sequence,data,sim,colors)

%% Plot orbit asteroids
% plot_orbits_asteroids(asteroid_names,colors)

%% delete the python file from this directory
delete('neo_api_function.py');