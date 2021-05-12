%% --------------------------------------------------------------------- %%
%% -------------------- EA Ast1 Ast2 Ast3 Ast4 Transfer ---------------- %%
%% --------------------------- ARCH 1+4, FLYBY ------------------------- %%
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

sim.case_name = 'ARCH ID 4: IMPULSIVE FLYBY ON EVERY ASTEROID WITH GA ON EARTH BEFORE GOING TO THE ASTEROIDS';

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

%% Spacecraft
sim.g0=9.81*(10^(-3)); %[km/s]
sim.Isp_mother=250;%s %Isp of mothercraft
sim.Isp_lander=230;%s Isp of landers
sim.dry_mass=40;%kg
sim.dry_mass_lander=4;%kg
% so that the starting dry mass is 40+4*4=56 kg

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
sim.soo_lim.TOF1_min = 100; % days
sim.soo_lim.TOF1_max = 3*365; % days
% TOF2 (3)
sim.soo_lim.TOF2_min = 50; % days
sim.soo_lim.TOF2_max = 3*365; % days
% TOF3 (4)
sim.soo_lim.TOF3_min = 50; % days
sim.soo_lim.TOF3_max = 3*365; % days
% TOF4 (5)
sim.soo_lim.TOF4_min = 50; % days
sim.soo_lim.TOF4_max = 3*365; % days
% Matrix of permutations (6)
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.soo_lim.permutations_low = 0.5; 
sim.soo_lim.permutations_up = data.HowMany + 0.4999;

% x = [MJD0,TOF1,TOF2,TOF3,TOF4,ID_permutation]
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed,  sim.soo_lim.TOF1_min,...
      sim.soo_lim.TOF2_min,sim.soo_lim.TOF3_min,...
      sim.soo_lim.TOF4_min,sim.soo_lim.permutations_low]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max,...
      sim.soo_lim.TOF2_max,sim.soo_lim.TOF3_max,...
      sim.soo_lim.TOF4_max,sim.soo_lim.permutations_up]; % Upper bound

% Constraint on C3 Launcher
sim.C3_max = 30; % km^2/s^2

%% Parallel pool
% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

%% Build the soo
FitnessFunction = @(x) ff_impulsive_soo_ARCH1plus4(x, data, sim); % Function handle to the fitness function
% FitnessFunction = @(x) ff_impulsive_soo_ARCH1plus4_mass(x, data, sim); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

%% Constraints
% sim.soo_constr.A = []; % linear inequality constraints
% sim.soo_constr.b = []; % linear inequality constraints
% sim.soo_constr.Aeq = []; % linear equality constraints
% sim.soo_constr.beq = []; % linear equality constraints
% sim.soo_constr.nonlcon = []; % linear equality constraints
% 
% % if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% % ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% sim.soo_constr.IntCon = [7];

% %% Options ga
% options = optimoptions(@ga);
% % options.PlotFcn = {@gaplotbestf};
% options.Display = 'iter';
% % A hybrid function is another minimization function that runs after the 
% % multiobjective genetic algorithm terminates
% % options.HybridFcn = @fmincon;
% 
% options.PopulationSize = 1000; 
% options.MaxGenerations = 200; 
% options.FunctionTolerance = 1e-6;
% options.MaxStallGenerations = 50;

%% Options ps
options = optimoptions('particleswarm');
options.HybridFcn = @fmincon;
options.SwarmSize = 1000; % Default is min(100,10*nvars),
options.MaxIterations = 300; %  Default is 200*nvars
options.MaxStallIterations = 50; % Default 20
options.Display = 'iter';
options.FunctionTolerance = 1e-6;
options.UseParallel = true;

%% Run Optimisation ga
% tic
% [x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,sim.soo_constr.A, ...
%     sim.soo_constr.b,sim.soo_constr.Aeq,sim.soo_constr.beq,sim.soo_bound.lb,sim.soo_bound.ub,...
%     sim.soo_constr.nonlcon,sim.soo_constr.IntCon,options);
% el_time_min_pp = toc/60;

%% Run Optimisation ps
numm=200;
for i = 1:numm
[x,Fval,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables...
    ,sim.soo_bound.lb,sim.soo_bound.ub,options);
el_time_min_pp = toc/60;

% Build solution structure
% set the knee as main solution
asteroid_sequence = data.PermutationMatrix(round(x(6)),:);
sol.ast_1 = asteroid_sequence(1);
sol.ast_2 = asteroid_sequence(2);
sol.ast_3 = asteroid_sequence(3);
sol.ast_4 = asteroid_sequence(4);
sol.MJD0 = x(1);
sol.dep_date = mjd20002date(sol.MJD0)';
sol.TOF_tot_D = x(2)+x(3)+x(4)+x(5);
sol.TOF_tot_Y = sol.TOF_tot_D/365;
sol.end_of_mission_date = mjd20002date(sol.MJD0+sol.TOF_tot_D)';
% sol.dV_tot = Fval(1); % case of dv optimisation without penalty
% sol.m_tot = Fval(1); % case of mass optimisation
sol.TOF1 = x(2);
sol.TOF2 = x(3);
sol.TOF3 = x(4);
sol.TOF4 = x(5);

sol = traj(sol,asteroid_sequence,data,sim);

AAA{i} = sol;
end

%% find min sol
for i=1:numm
    dv_vect(i)=AAA{i}.dV_tot;
end
idx_min=find(min(dv_vect)==dv_vect);
sol=AAA{idx_min};

%% Mass Consumption for High Thrust Impulsive Case
% g0 = 9.81; %m/s^2
% % https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
% Isp = 230; %s 
% m_dry = 100; %kg
% m_prop = m_dry*(exp(sol.dV_tot*1e3/(g0*Isp)) - 1); %kg
%% Plot trajectories
% sol = plot_mission_4neo_flyby_GA(sol,asteroid_sequence,data,sim,colors)

sol = plot_mission_4neo_flyby(sol,asteroid_sequence,data,sim,colors)

%% Plot Angles
figure()
plot(sol.SCtime,sol.angles.SCA)
xlabel('time [s]'); ylabel('Sun Conjunction Angle [rad]');
hold on
plot(sol.SCtime,sol.angles.SolarConjunction)
xlabel('time [s]'); ylabel('Sun Conjunction Angle [rad]');
figure()
plot(sol.SCtime,sol.angles.SAA)
xlabel('time [s]'); ylabel('Sun Aspect Angle [rad]');
figure()
plot(sol.SCtime,sol.angles.EVA)
xlabel('time [s]'); ylabel('Earth Visibility Angle [rad]');
%% Plot orbit asteroids
% plot_orbits_asteroids(asteroid_names,colors)

%% delete the python file from this directory
delete('neo_api_function.py');