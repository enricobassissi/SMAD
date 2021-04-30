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

sim.case_name = 'ARCH ID 1: IMPULSIVE FLYBY ON EACH ASTEROID';
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

% Number of possible combination of 4 asteroids among the ones in the list
data.HowMany = factorial(length(data.asteroid_names)) / factorial(length(data.asteroid_names) - 4);
[data.PermutationMatrix, ~] = permnUnique(data.asteroid_names, 4);

%% uNEO
% try 
%     load('data.mat')
% catch
    % if the asteroid have changed, run the find_eph_neo below, it takes about 1 min
    [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);
%     save('data.mat', data);
% end
%% Boundaries
% Departure dates (1)
sim.bound.date_ed = [2022, 1, 1, 0, 0, 0];
sim.bound.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.bound.mjd2000_ed = date2mjd2000(sim.bound.date_ed);
sim.bound.mjd2000_ld = date2mjd2000(sim.bound.date_ld);
% TOF1 (2)
sim.bound.TOF1_min = 200; % days
sim.bound.TOF1_max = 3*365; % days
% TOF2 (3)
sim.bound.TOF2_min = 50; % days
sim.bound.TOF2_max = 3*365; % days
% TOF3 (4)
sim.bound.TOF3_min = 50; % days
sim.bound.TOF3_max = 3*365; % days
% TOF4 (5)
sim.bound.TOF4_min = 50; % days
sim.bound.TOF4_max = 3*365; % days
% Matrix of permutations (6)
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.bound.permutations_low = 0.5; 
sim.bound.permutations_up = data.HowMany + 0.4999;

% x = [MJD0,TOF1,TOF2,TOF3,TOF4,ID_permutation,]
sim.bound.lb = [sim.bound.mjd2000_ed, sim.bound.TOF1_min,...
    sim.bound.TOF2_min,sim.bound.TOF3_min,sim.bound.TOF4_min,...
    sim.bound.permutations_low]; % Lower bound
sim.bound.ub = [sim.bound.mjd2000_ld, sim.bound.TOF1_max,...
    sim.bound.TOF2_max,sim.bound.TOF3_max,...
    sim.bound.TOF4_max,sim.bound.permutations_up]; % Upper bound

% Constraint on C3 Launcher
sim.C3_max = 20; % km^2/s^2

%% Options

% Parallel pool
% Open the parallel pool
% par_pool = gcp; 
% if isempty(par_pool)
%     poolsize = 0;
% else
%     poolsize = par_pool.NumWorkers;
% end
% 
% options.UseParallel = true;
% 
% options.UseParallel = false;

%% Build the moo
FitnessFunction = @(x) ff_impulsive_moo_ps_dVtot_dVrel(x, data, sim); % Function handle to the fitness function
numberOfVariables = length(sim.bound.ub); % Number of decision variables

%% Mopso parameter
params.Np = 1000;        % Population size
params.Nr = 200;        % Repository size
params.maxgen = 200;    % Maximum number of generations
params.W = 0.4;         % Inertia weight
params.C1 = 2;          % Individual confidence factor
params.C2 = 2;          % Swarm confidence factor
params.ngrid = 20;      % Number of grids in each dimension
params.maxvel = 5;      % Maxmium vel in percentage
params.u_mut = 0.5;     % Uniform mutation percentage


MultiObj.fun = FitnessFunction;
MultiObj.nVar = length(sim.bound.ub);
MultiObj.var_min = sim.bound.lb;
MultiObj.var_max = sim.bound.ub;
MultiObj.Obj1 = '$\Delta V_{tot}$ [km/s]';
MultiObj.Obj2 = '$\Delta V_{rel \ passage}$ [km/s]';

%% MOPSO
REP = MOPSO(params,MultiObj);

%% Find the knee solution
[knee_idx, d] = find_knee_solution(REP.pos_fit);

% Plot Pareto Plot
figure('Name','PS MO Pareto Plot')
title('Pareto Points in Parameter Space')
h_pp = plot(REP.pos_fit(:,1),REP.pos_fit(:,2),'o','Color',colors(1,:));
hold on
h_kpp = plot(REP.pos_fit(knee_idx,1),REP.pos_fit(knee_idx,2),'o','Color',colors(2,:));
xlabel(MultiObj.Obj1)
ylabel(MultiObj.Obj2)
legend([h_pp,h_kpp],'Sub-Optim Sol','Knee Sol')
clearvars h_pp h_kpp

%% find min deltav sol
knee_idx = find(min(REP.pos_fit(:,1))==REP.pos_fit(:,1));

%% Build solution structure
% set the knee as main solution
asteroid_sequence = data.PermutationMatrix(round(REP.pos(knee_idx,6)),:);
sol.ast_1 = asteroid_sequence(1);
sol.ast_2 = asteroid_sequence(2);
sol.ast_3 = asteroid_sequence(3);
sol.ast_4 = asteroid_sequence(4);
sol.TOF1 = REP.pos(knee_idx,2);
sol.TOF2 = REP.pos(knee_idx,3);
sol.TOF3 = REP.pos(knee_idx,4);
sol.TOF4 = REP.pos(knee_idx,5);
sol.TOF_tot = sol.TOF1+sol.TOF2+sol.TOF3+sol.TOF4;
sol.MJD0 = REP.pos(knee_idx,1);
sol.dep_date = mjd20002date(sol.MJD0)';
sol.end_of_mission_date = mjd20002date(sol.MJD0+sol.TOF_tot)';
sol.dV_tot = REP.pos_fit(knee_idx,1);
sol.dV_rel = REP.pos_fit(knee_idx,2);

%% Mass Consumption for High Thrust Impulsive Case
g0 = 9.81; %m/s^2
% https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
Isp = 230; %s 
m_dry = 100; %kg
m_prop = m_dry*(exp(sol.dV_tot*1e3/(g0*Isp)) - 1); %kg
%% Plot trajectories
sol = plot_mission_4neo_flyby(sol,asteroid_sequence,data,sim,colors)

%% Plot orbit asteroids
% plot_orbits_asteroids(asteroid_names,colors)

%% delete the python file from this directory
delete('neo_api_function.py');