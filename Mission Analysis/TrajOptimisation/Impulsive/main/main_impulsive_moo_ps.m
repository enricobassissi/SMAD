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

%% clear and close
clear; close all; clc;
% addpath time
% addpath function
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
%% Initializing the Environment
AU = astroConstants(2);
muSun = astroConstants(4);
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

% %% INTRO ADIMENSIONALISATION
% sim.mu = 1.32712440017987e11; % Sun planetary constant (mu = mass * G) (from DE405) [km^3/s^2]
% sim.DU = 149597870.691; % Distance Unit = Astronomical Unit (AU) (from DE405) [km]
% sim.TU = (sim.DU^3/sim.mu)^0.5; % Time Unit
% sim.mu = 1;

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
% Departure dates
sim.bound.date_ed = [2022, 1, 1, 0, 0, 0];
sim.bound.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.bound.mjd2000_ed = date2mjd2000(sim.bound.date_ed);
sim.bound.mjd2000_ld = date2mjd2000(sim.bound.date_ld);
% TOF1
sim.bound.TOF1_min = 200; % days
sim.bound.TOF1_max = 3*365; % days
% % Launcher velocity given and angles
% sim.bound.v_inf_magn_min = 0;
% sim.bound.v_inf_magn_max = sqrt(40); % c3 = 40 km/s^2
% sim.bound.alpha_min = deg2rad(0);
% sim.bound.alpha_max = deg2rad(360);
% sim.bound.beta_min = deg2rad(0);
% sim.bound.beta_max = deg2rad(360);
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
sim.bound.permutations_up = data.HowMany + 0.4999;
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
sim.bound.lb = [sim.bound.mjd2000_ed, sim.bound.TOF1_min , sim.bound.bt1_min,...
      sim.bound.TOF2_min,sim.bound.permutations_low,sim.bound.bt2_min,...
      sim.bound.TOF3_min,sim.bound.bt3_min,sim.bound.TOF4_min]; % Lower bound
sim.bound.ub = [sim.bound.mjd2000_ld, sim.bound.TOF1_max, sim.bound.bt1_max,...
      sim.bound.TOF2_max,sim.bound.permutations_up,sim.bound.bt2_max,...
      sim.bound.TOF3_max,sim.bound.bt3_max,sim.bound.TOF4_max]; % Upper bound
 % Constraint on C3 Launcher
sim.C3_max = 40; % km^2/s^2
%% Constraints
sim.constr.A = []; % linear inequality constraints
sim.constr.b = []; % linear inequality constraints
sim.constr.Aeq = []; % linear equality constraints
sim.constr.beq = []; % linear equality constraints
sim.constr.nonlcon = []; % linear equality constraints

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
FitnessFunction = @(x) ff_impulsive_moo_ps(x, data, sim); % Function handle to the fitness function
numberOfVariables = length(sim.bound.ub); % Number of decision variables


%% Mopso parameter
params.Np = 200;        % Population size
params.Nr = 200;        % Repository size
params.maxgen = 100;    % Maximum number of generations
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
%% MOPSO
REP = MOPSO(params,MultiObj);
