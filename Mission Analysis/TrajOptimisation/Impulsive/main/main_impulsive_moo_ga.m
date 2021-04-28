%% -------------------------------------------------------- %%
%% ----------- EA Ast1 Ast2 Ast3 Ast4 Transfer ------------ %%
%% -------------------------------------------------------- %%
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

%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\Impulsive\main');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'main');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% INTRO ADIMENSIONALISATION
AU = astroConstants(2);
muSun = astroConstants(4);

sim.mu = 1.32712440017987e11; % Sun planetary constant (mu = mass * G) (from DE405) [km^3/s^2]
sim.DU = 149597870.691; % Distance Unit = Astronomical Unit (AU) (from DE405) [km]
sim.TU = (sim.DU^3/sim.mu)^0.5; % Time Unit
sim.mu = 1;

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%%
% Asteroids
% data extraction section
asteroid_names = ["2006HX57";"2008XU2";"2008KN11";"2012SY49";"2012QD8";"2020UE";...
                  "2006SC";"2005WG57";"2012BY1"];

% Number of possible combination of 4 asteroids among the ones in the list
HowMany = factorial(length(asteroid_names)) / factorial(length(asteroid_names) - 4);
[PermutationMatrix, ~] = permnUnique(asteroid_names, 4);

%% Boundaries
% Departure dates
sim.bound.date_ed = [2022, 1, 1, 0, 0, 0];
sim.bound.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.bound.mjd2000_ed = date2mjd2000(sim.bound.date_ed);
sim.bound.mjd2000_ld = date2mjd2000(sim.bound.date_ld);
% TOF1
sim.bound.TOF1_min = 200; % days
sim.bound.TOF1_max = 3*365; % days
% Launcher velocity given and angles
sim.bound.v_inf_magn_min = 0;
sim.bound.v_inf_magn_max = sqrt(40); % c3 = 40 km/s^2
sim.bound.alpha_min = deg2rad(0);
sim.bound.alpha_max = deg2rad(360);
sim.bound.beta_min = deg2rad(0);
sim.bound.beta_max = deg2rad(360);
% Buffer time 1
sim.bound.bt1_min = 30;
sim.bound.bt1_max = 180;
% TOF2
sim.bound.TOF2_min = 50; % days
sim.bound.TOF2_max = 3*365; % days
% Matrix of permutations
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.bound.permutations_low = 0.5; 
sim.bound.permutations_up = HowMany + 0.4999;
% Buffer time 2
sim.bound.bt2_min = 30;
sim.bound.bt2_max = 180;
% TOF3
sim.bound.TOF3_min = 50; % days
sim.bound.TOF3_max = 3*365; % days
% Buffer time 3 
sim.bound.bt3_min = 30;
sim.bound.bt3_max = 180;
% TOF4
sim.bound.TOF4_min = 50; % days
sim.bound.TOF4_max = 3*365; % days

% x = [MJD0,TOF1,v_inf_magn,aplha,beta,buffer_time,TOF2,ID_permutation,...
%      buffer_time2,TOF3,buffer_time3,TOF4]
sim.bound.lb = [sim.bound.mjd2000_ed, sim.bound.TOF1_min, sim.bound.v_inf_magn_min,...
      sim.bound.alpha_min, sim.bound.beta_min, sim.bound.bt1_min,...
      sim.bound.TOF2_min,sim.bound.permutations_low,sim.bound.bt2_min,...
      sim.bound.TOF3_min,sim.bound.bt3_min,sim.bound.TOF4_min]; % Lower bound
sim.bound.ub = [sim.bound.mjd2000_ld, sim.bound.TOF1_max, sim.bound.v_inf_magn_max,...
      sim.bound.alpha_max, sim.bound.beta_max, sim.bound.bt1_max,...
      sim.bound.TOF2_max,sim.bound.permutations_up,sim.bound.bt2_max,...
      sim.bound.TOF3_max,sim.bound.bt3_max,sim.bound.TOF4_max]; % Upper bound
  
%% Constraints
sim.constr.A = []; % linear inequality constraints
sim.constr.b = []; % linear inequality constraints
sim.constr.Aeq = []; % linear equality constraints
sim.constr.beq = []; % linear equality constraints
sim.constr.nonlcon = []; % linear equality constraints

%% Options
options = optimoptions(@gamultiobj);
% options.PlotFcn = @gaplotpareto;
options.Display = 'iter';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
% options.HybridFcn = 'fgoalattain';

options.PopulationSize = 3; % ideal 1000
options.ParetoFraction = 0.7;
options.MaxGenerations = 3; % ideal 100
options.FunctionTolerance = 1e-6;
options.MaxStallGenerations = 3;

% Parallel pool
% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;

options.UseParallel = false;

%% Build the moo
FitnessFunction = @(x) ff_neo_perm(x, PermutationMatrix); % Function handle to the fitness function
numberOfVariables = length(sim.bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = gamultiobj(FitnessFunction,numberOfVariables,sim.constr.A, ...
    sim.constr.b,sim.constr.Aeq,sim.constr.beq,sim.bound.lb,sim.bound.ub,...
    sim.constr.nonlcon,options);
el_time_min_pp = toc/60;

%% Find the knee solution
[knee_idx, d] = find_knee_solution(Fval);

% Plot Pareto Plot
figure('Name','GA MO Pareto Plot')
title('Pareto Points in Parameter Space')
h_pp = plot(Fval(:,1),Fval(:,2),'o','Color',colors(1,:));
hold on
h_kpp = plot(Fval(knee_idx,1),Fval(knee_idx,2),'o','Color',colors(2,:));
xlabel('$Obj_1: \ \Delta V$ [km/s]')
ylabel('$Obj_2: \ TOF$ [d]')
legend([h_pp,h_kpp],'Sub-Optim Sol','Knee Sol')
clearvars h_pp h_kpp

%% Build solution structure
% set the knee as main solution
asteroid_sequence = PermutationMatrix(round(x(knee_idx,8)),:);
sol.ast_1 = asteroid_sequence(1);
sol.ast_2 = asteroid_sequence(2);
sol.ast_3 = asteroid_sequence(3);
sol.ast_4 = asteroid_sequence(4);
sol.MJD0 = x(knee_idx,1);
sol.dep_date = mjd20002date(sol.MJD0)';
sol.end_of_mission_date = mjd20002date(sol.MJD0+Fval(knee_idx,2))';
sol.dV_tot = Fval(knee_idx,1);
sol.TOF_tot = Fval(knee_idx,2);
sol.TOF1 = x(knee_idx,2);
sol.buffer_time1 = x(knee_idx,6);
sol.TOF2 = x(knee_idx,7);
sol.buffer_time2 = x(knee_idx,9);
sol.TOF3 = x(knee_idx,10);
sol.buffer_time3 = x(knee_idx,11);
sol.TOF4 = x(knee_idx,12);
sol.v_inf_magn = x(knee_idx,3);
sol.v_inf_alpha = rad2deg(x(knee_idx,4));
sol.v_inf_beta = rad2deg(x(knee_idx,5));

%% Mass Consumption for High Thrust Impulsive Case
g0 = 9.81; %m/s^2
% https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
Isp = 230; %s 
m_dry = 100; %kg
m_prop = m_dry*(exp(sol.dV_tot*1e3/(g0*Isp)) - 1); %kg
%% Plot trajectories
sol = plot_mission_4neo(sol,colors,asteroid_sequence)

%% Plot orbit asteroids
% plot_orbits_asteroids(asteroid_names,colors)

%% delete the python file from this directory
delete('neo_api_function.py');