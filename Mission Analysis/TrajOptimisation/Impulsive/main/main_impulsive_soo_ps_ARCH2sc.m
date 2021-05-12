%% -------------------------------------------------------- %%
%% ----------- EA Ast1 Ast2 Ast3 Ast4 Transfer ------------ %%
%% ------------------- ARCH ID 3: 2 SC -------------------- %%
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

sim.case_name = 'ARCH ID 3: 2 SC, EACH IMPULSIVE FLYBY ON 2 ASTEROIDS';
% In this architecture you have all the permutation of 2 asteorids among
% the list of 9
% Then the 2nd spacecraft can combine freely the other targets, except from
% the 2 asteroids that have already been set for the 1st spacecraft

% %% INTRO ADIMENSIONALISATION
% sim.mu = 1.32712440017987e11; % Sun planetary constant (mu = mass * G) (from DE405) [km^3/s^2]
% sim.DU = 149597870.691; % Distance Unit = Astronomical Unit (AU) (from DE405) [km]
% sim.TU = (sim.DU^3/sim.mu)^0.5; % Time Unit
% sim.mu = 1;

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

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%% Asteroids
AU = astroConstants(2);
muSun = astroConstants(4);

% data extraction section
data.asteroid_names = ["2006HX57";"2008XU2";"2008KN11";"2012SY49";"2012QD8";"2020UE";...
                  "2006SC";"2005WG57";"2012BY1"];

% Number of possible combination of 2 asteroids among the ones in the list
data.HowMany = factorial(length(data.asteroid_names)) / factorial(length(data.asteroid_names) - 2);
[data.PermutationMatrix, ~] = permnUnique(data.asteroid_names, 2);

% the left asteroids permutations
data.HowMany2 = factorial(length(data.asteroid_names)-2) / factorial((length(data.asteroid_names)-2) - 2);
%% uNEO
% try 
clearvars data
load('data_2SC.mat');
% catch
    % if the asteroid have changed, run the find_eph_neo below, it takes about 1 min
%     [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);
%     save('data.mat', data);
% end

%% Boundaries
% Nomenclature
% quantities for SpaceCraft1 will have numbers (1,2,...)
% quantities for SC 2 will have letters (a,b,...)

% Departure dates (1), departure time for both the sc
sim.soo_lim.date_ed = [2022, 1, 1, 0, 0, 0];
sim.soo_lim.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.soo_lim.mjd2000_ed = date2mjd2000(sim.soo_lim.date_ed);
sim.soo_lim.mjd2000_ld = date2mjd2000(sim.soo_lim.date_ld);
% TOF1 (2), tof sc1 to 1st asteroid
sim.soo_lim.TOF1_min = 200; % days
sim.soo_lim.TOF1_max = 3*365; % days
% TOF2 (3), tof sc1 to 2nd asteroid
sim.soo_lim.TOF2_min = 50; % days
sim.soo_lim.TOF2_max = 3*365; % days
% TOFa (4), tof sc2 to 1st asteroid
sim.soo_lim.TOF3_min = 200; % days
sim.soo_lim.TOF3_max = 3*365; % days
% TOFb (5), tof sc2 to 2nd asteroid
sim.soo_lim.TOF4_min = 50; % days
sim.soo_lim.TOF4_max = 3*365; % days
% Matrix of permutations 1 (6)
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.soo_lim.permutations_low = 0.5; 
sim.soo_lim.permutations_up = data.HowMany + 0.4999;
% Matrix of permutations 2 (7)
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.soo_lim.permutations2_low = 0.5; 
sim.soo_lim.permutations2_up = data.HowMany2 + 0.4999;

% x = [MJD0,TOF1,TOF2,TOF3,TOF4,ID_permutation1,ID_permutation2]
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed, sim.soo_lim.TOF1_min,...
      sim.soo_lim.TOF2_min,sim.soo_lim.TOF3_min,sim.soo_lim.TOF4_min,...
      sim.soo_lim.permutations_low,sim.soo_lim.permutations2_low]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max,...
      sim.soo_lim.TOF2_max,sim.soo_lim.TOF3_max,sim.soo_lim.TOF4_max,...
      sim.soo_lim.permutations_up,sim.soo_lim.permutations2_up]; % Upper bound

% Constraint on C3 Launcher
sim.C3_max = 30; % km^2/s^2

%% Options
options = optimoptions('particleswarm');
options.HybridFcn = @fmincon;
options.SwarmSize = 1000; % Default is min(100,10*nvars),
options.MaxIterations = 200; %  Default is 200*nvars
options.MaxStallIterations = 70; % Default 20
options.Display = 'iter';
options.FunctionTolerance = 1e-6;

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
FitnessFunction = @(x) ff_impulsive_soo_ARCH2sc(x, data, sim); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

%% Run the soo
tic
[x,Fval,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables...
    ,sim.soo_bound.lb,sim.soo_bound.ub,options);
el_time_min_pp = toc/60;

%% Build solution structure
sol.MJD0 = x(1);
sol.dep_date = mjd20002date(sol.MJD0)';
% 1st sc stuff
asteroid_sequence = data.PermutationMatrix(round(x(6)),:);
sol.ast_1 = asteroid_sequence(1);
sol.ast_2 = asteroid_sequence(2);
sol.TOF1 = x(2);
sol.TOF2 = x(3);
sol.TOF_tot_D_sc1 = sol.TOF1+sol.TOF2;
sol.TOF_tot_Y_sc1 = sol.TOF_tot_D_sc1/365;
sol.end_of_mission_date_sc1 = mjd20002date(sol.MJD0+sol.TOF_tot_D_sc1)';
% sol.dV_tot = Fval(1);

% 2nd sc stuff
TF = contains(data.asteroid_names,asteroid_sequence);
data.not_asteroid_sequence = data.asteroid_names(~TF);
clearvars TF
% HowMany_for2ndSC = factorial(length(not_asteroid_sequence)) / factorial(length(not_asteroid_sequence) - 2);
[data.PermutationMatrix_SC2, ~] = permnUnique(data.not_asteroid_sequence, 2);
sol.ast_a = data.PermutationMatrix_SC2(round(x(7)),1);
sol.ast_b = data.PermutationMatrix_SC2(round(x(7)),2);
sol.TOFa = x(4);
sol.TOFb = x(5);
sol.TOF_tot_D_sc2 = sol.TOFa+sol.TOFb;
sol.TOF_tot_Y_sc2 = sol.TOF_tot_D_sc2/365;
sol.end_of_mission_date_sc2 = mjd20002date(sol.MJD0+sol.TOF_tot_D_sc2)';

sol.end_of_mission_date_overall = max(sol.end_of_mission_date_sc1,sol.end_of_mission_date_sc2);

% %% Mass Consumption for High Thrust Impulsive Case
% g0 = 9.81; %m/s^2
% % https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
% Isp = 230; %s 
% m_dry = 100; %kg
% m_prop = m_dry*(exp(sol.dV_tot*1e3/(g0*Isp)) - 1); %kg

%% Plot trajectories
sol = plot_mission_4neo_flyby_ARCH2sc(sol,data,sim,colors)
[sol_dates] = sol_to_dates_of_mission(sol,'2SC_FB')

%% Plot orbit asteroids
% plot_orbits_asteroids(asteroid_names,colors)

%% delete the python file from this directory
delete('neo_api_function.py');