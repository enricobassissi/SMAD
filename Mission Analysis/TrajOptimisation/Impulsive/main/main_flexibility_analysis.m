%% -------------------------------------------------------- %%
%% ----------- EA Ast1 Ast2 Ast3 Ast4 Transfer ------------ %%
%% ---------------- FLEXIBILITY ANALYSIS ------------------ %%
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
          0    174  157;... % (12) PURE TEAL
          0    0    0]./255; % (13) BLACK

sim.case_name = 'ARCH ID 2: IMPULSIVE FLYBY ON EVERY ASTEROID';

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
%     if the asteroid have changed, run the find_eph_neo below, it takes about 1 min
[data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

%% delete the python file from this directory
delete('neo_api_function.py');

%% Boundaries
% Departure dates (1)
sim.soo_lim.date_ed = [2022, 1, 1, 0, 0, 0];
sim.soo_lim.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.soo_lim.mjd2000_ed = date2mjd2000(sim.soo_lim.date_ed);
sim.soo_lim.mjd2000_ld = date2mjd2000(sim.soo_lim.date_ld);
% TOF1 (2)
sim.soo_lim.TOF1_min = 200; % days
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
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed, sim.soo_lim.TOF1_min,...
      sim.soo_lim.TOF2_min,sim.soo_lim.TOF3_min,...
      sim.soo_lim.TOF4_min,sim.soo_lim.permutations_low]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max,...
      sim.soo_lim.TOF2_max,sim.soo_lim.TOF3_max,...
      sim.soo_lim.TOF4_max,sim.soo_lim.permutations_up]; % Upper bound

% Constraint on C3 Launcher
sim.C3_max = 20; % km^2/s^2

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
FitnessFunction = @(x) ff_impulsive_soo_ARCH1plus4(x, data, sim); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables...
    ,sim.soo_bound.lb,sim.soo_bound.ub,options);
el_time_min_pp = toc/60;

%% Build solution structure
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
sol.dV_tot = Fval(1);
sol.TOF1 = x(2);
sol.TOF2 = x(3);
sol.TOF3 = x(4);
sol.TOF4 = x(5);

%% Mass Consumption for High Thrust Impulsive Case
g0 = 9.81; %m/s^2
% https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
Isp = 230; %s 
m_dry = 100; %kg
m_prop = m_dry*(exp(sol.dV_tot*1e3/(g0*Isp)) - 1); %kg
%% Plot trajectories
sol = plot_mission_4neo_flyby(sol,asteroid_sequence,data,sim,colors)

%% ------------------------------------------------ %%
%% ----------- FLEXIBILITY ANALYSIS --------------- %%
%% ------------ OVER LAUNCH DATE ------------------ %%
%% ------------------------------------------------ %%
%% Synodic period Earth - Asteroid 1
% it loose of sense if the asteroid sequence changes
% [kep_EA,ksun] = uplanet(sol.MJD0, 3);
% [kep_ast_1] = uNEO2(sol.MJD0+x(2),sol.ast_1,data);
% syn_periods_EA_ast1 = synodic_period(kep_EA(1),kep_ast_1(1))/(3600*24); % days

% let's make arbitrary that the max delay time is the one of Rosetta
% mission, one of the most delayed mission
delay_time = 730;

%% launch windows vector
% that now it becomes a fixed variable and not anymore a degree of freedom
% if i have problem for the launch, if i posticipate the launch is it a
% problem? meaning how does the costs (dV, TOF) changes?
number_of_points = 10; % how many points do you want? fine or coarse grid
lw_vector = linspace(sol.MJD0,sol.MJD0+delay_time,delay_time/number_of_points);

%% Permutations among the 4 asteroid fixed
% Number of possible combination of 4 asteroids among the ones in the list
data.HowMany_flex = factorial(length(asteroid_sequence)) / factorial(length(asteroid_sequence) - 4);
[data.PermutationMatrix_flex, ~] = permnUnique(asteroid_sequence, 4);

%% Boundaries
% the sequence now is the one dictated by the best sol found
% TOF1 (1)
sim.soo_lim.TOF1_min = 200; % days
sim.soo_lim.TOF1_max = 3*365; % days
% TOF2 (2)
sim.soo_lim.TOF2_min = 50; % days
sim.soo_lim.TOF2_max = 3*365; % days
% TOF3 (3)
sim.soo_lim.TOF3_min = 50; % days
sim.soo_lim.TOF3_max = 3*365; % days
% TOF4 (4)
sim.soo_lim.TOF4_min = 50; % days
sim.soo_lim.TOF4_max = 3*365; % days
% Matrix of permutations (5)
% to use round in the code... so we have same probility to be rounded to
% the first or to the last element in the matrix as in the middle elements!
sim.soo_lim.permutations_low = 0.5; 
sim.soo_lim.permutations_up = data.HowMany_flex + 0.4999;

% x = [TOF1,TOF2,TOF3,TOF4,ID_Permutation]
sim.soo_bound.lb = [sim.soo_lim.TOF1_min,sim.soo_lim.TOF2_min,...
    sim.soo_lim.TOF3_min,sim.soo_lim.TOF4_min,sim.soo_lim.permutations_low]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.TOF1_max,sim.soo_lim.TOF2_max,...
    sim.soo_lim.TOF3_max,sim.soo_lim.TOF4_max,sim.soo_lim.permutations_up]; % Upper bound

% Constraint on C3 Launcher
sim.C3_max = 20; % km^2/s^2

%% PS Options
options = optimoptions('particleswarm');
% options.HybridFcn = @fmincon;
options.SwarmSize = 500; % Default is min(100,10*nvars),
options.MaxIterations = 100; %  Default is 200*nvars
options.MaxStallIterations = 50; % Default 20
options.Display = 'final';
options.FunctionTolerance = 1e-6;

% Parallel pool - Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;

%% SOO but looping over the dep date
% and saving each result on dV_tot and TOF_tot
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

% initialize
dV_flex = zeros(length(lw_vector),1);
TOF_flex = zeros(length(lw_vector),1);
% loop over MJD0 for sensitivity analysis over departure date
for i = 1:length(lw_vector)
    sim.MJD0_flexibility = lw_vector(i);
    FitnessFunction = @(x) ff_impulsive_flexibility(x, data, sim); % Function handle to the fitness function
    [x_flex,Fval_flex,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables,...
        sim.soo_bound.lb,sim.soo_bound.ub,options);
    
    % save the result
    flex{i}.asteroid_sequence = data.PermutationMatrix_flex(round(x_flex(5)),:);
    flex{i}.ast_1 = flex{i}.asteroid_sequence(1);
    flex{i}.ast_2 = flex{i}.asteroid_sequence(2);
    flex{i}.ast_3 = flex{i}.asteroid_sequence(3);
    flex{i}.ast_4 = flex{i}.asteroid_sequence(4);
    flex{i}.MJD0 = sim.MJD0_flexibility;
    flex{i}.dep_date = mjd20002date(sim.MJD0_flexibility)';
    flex{i}.TOF_tot_D = x_flex(1)+x_flex(2)+x_flex(3)+x_flex(4);
    flex{i}.TOF_tot_Y = flex{i}.TOF_tot_D/365;
    flex{i}.end_of_mission_date = mjd20002date(sim.MJD0_flexibility+flex{i}.TOF_tot_D)';
    flex{i}.dV_tot = Fval_flex;
    flex{i}.TOF1 = x_flex(1);
    flex{i}.TOF2 = x_flex(2);
    flex{i}.TOF3 = x_flex(3);
    flex{i}.TOF4 = x_flex(4);

end
clearvars i

%% plot the result
for i=1:length(lw_vector)
    dV_flex(i,1) = flex{i}.dV_tot;
    TOF_flex(i,1) = flex{i}.TOF_tot_D;
end
clearvars i

figure('Name','Flexibility over Departure Time and Asteroids Sequence Order')
yyaxis left
plot(lw_vector',dV_flex)
ylabel('$\Delta V$ [km/s]')
hold on
yyaxis right
plot(lw_vector',TOF_flex)
ylabel('TOF [d]')
xlabel('Departure MJD2000')
% plot(sol.MJD0,sol.dV_tot,'o')

%% ------------------------------------------------ %%
%% ----------- FLEXIBILITY ANALYSIS --------------- %%
%% ---------- Let's change one asteroid ----------- %%
%% ------------------------------------------------ %%
%% Permutations among the 4 asteroid fixed
% Number of possible combination of 4 asteroids among the ones in the list
data.HowMany_flex = factorial(length(asteroid_sequence)) / factorial(length(asteroid_sequence) - 4);
[data.PermutationMatrix_flex, ~] = permnUnique(asteroid_sequence, 4);

TF = contains(data.asteroid_names,asteroid_sequence);
not_asteroid_sequence = data.asteroid_names(~TF);
clearvars TF

data.NumberOfOtherAsteroids = length(not_asteroid_sequence);

for i=1:data.NumberOfOtherAsteroids % rows
    for j=1:length(asteroid_sequence) % columns
        new_asteroid_sequence = asteroid_sequence;
        new_asteroid_sequence(j) = not_asteroid_sequence(i);
        data.AstFlexSequenceMatrix{i,j} = permnUnique(new_asteroid_sequence, 4);
    end
end
clearvars i j

%% Boundaries
% Departure dates (1)
sim.soo_lim.date_ed = [2022, 1, 1, 0, 0, 0];
sim.soo_lim.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.soo_lim.mjd2000_ed = date2mjd2000(sim.soo_lim.date_ed);
sim.soo_lim.mjd2000_ld = date2mjd2000(sim.soo_lim.date_ld);
% TOF1 (2)
sim.soo_lim.TOF1_min = 200; % days
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
sim.soo_lim.permutations_up = data.HowMany_flex + 0.4999;

% x = [MJD0,TOF1,TOF2,TOF3,TOF4,ID_permutation]
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed, sim.soo_lim.TOF1_min,...
      sim.soo_lim.TOF2_min,sim.soo_lim.TOF3_min,...
      sim.soo_lim.TOF4_min,sim.soo_lim.permutations_low]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max,...
      sim.soo_lim.TOF2_max,sim.soo_lim.TOF3_max,...
      sim.soo_lim.TOF4_max,sim.soo_lim.permutations_up]; % Upper bound
% Constraint on C3 Launcher
sim.C3_max = 20; % km^2/s^2

%% PS Options
options = optimoptions('particleswarm');
% options.HybridFcn = @fmincon;
options.SwarmSize = 500; % Default is min(100,10*nvars),
options.MaxIterations = 100; %  Default is 200*nvars
options.MaxStallIterations = 50; % Default 20
options.Display = 'final';
options.FunctionTolerance = 1e-6;

% Parallel pool - Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;

%% SOO but looping over the dep date
% and saving each result on dV_tot and TOF_tot
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

% initialize
idx = 0;
% dV_flex = zeros(length(lw_vector),1);
% TOF_flex = zeros(length(lw_vector),1);
% loop over MJD0 for sensitivity analysis over departure date
for i=1:data.NumberOfOtherAsteroids % rows
    for j=1:length(asteroid_sequence) % columns
        data.PermutCycle = data.AstFlexSequenceMatrix{i,j};
        FitnessFunction = @(x) ff_impulsive_flexibility_ast_change(x, data, sim); % Function handle to the fitness function
        [x_flex,Fval_flex,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables,...
            sim.soo_bound.lb,sim.soo_bound.ub,options);

        % save the result
        flex_order{i,j}.asteroid_sequence = data.PermutCycle(round(x_flex(6)),:);
        flex_order{i,j}.ast_1 = flex_order{i,j}.asteroid_sequence(1);
        flex_order{i,j}.ast_2 = flex_order{i,j}.asteroid_sequence(2);
        flex_order{i,j}.ast_3 = flex_order{i,j}.asteroid_sequence(3);
        flex_order{i,j}.ast_4 = flex_order{i,j}.asteroid_sequence(4);
        flex_order{i,j}.MJD0 = x_flex(1);
        flex_order{i,j}.dep_date = mjd20002date(x_flex(1))';
        flex_order{i,j}.TOF_tot_D = x_flex(2)+x_flex(3)+x_flex(4)+x_flex(5);
        flex_order{i,j}.TOF_tot_Y = flex_order{i,j}.TOF_tot_D/365;
        flex_order{i,j}.end_of_mission_date = mjd20002date(x_flex(1)+flex_order{i,j}.TOF_tot_D)';
        flex_order{i,j}.dV_tot = Fval_flex;
        flex_order{i,j}.TOF1 = x_flex(2);
        flex_order{i,j}.TOF2 = x_flex(3);
        flex_order{i,j}.TOF3 = x_flex(4);
        flex_order{i,j}.TOF4 = x_flex(5);
        
        idx = idx+1;
        fprintf('Cycle %d \n',idx);
    end
end
clearvars i j idx

%% plot the result
for i=1:data.NumberOfOtherAsteroids % rows
    for j=1:length(asteroid_sequence) % columns
        dV_flex_order(i,j) = flex_order{i,j}.dV_tot;
        TOF_flex_order(i,j) = flex_order{i,j}.TOF_tot_D;
    end
end
clearvars i j

figure('Name','Changing and optimising the asteroid sequence order')
scatter(dV_flex_order(:,1),TOF_flex_order(:,1),[],colors(1,:),'filled',...
    'DisplayName','Ast 1 Variations and Reordering'); % all changes on first ast
hold on
scatter(dV_flex_order(:,2),TOF_flex_order(:,2),[],colors(2,:),'filled',...
    'DisplayName','Ast 2 Variations and Reordering'); % all changes on 2nd ast
scatter(dV_flex_order(:,3),TOF_flex_order(:,3),[],colors(4,:),'filled',...
    'DisplayName','Ast 3 Variations and Reordering'); % all changes on 3rd ast
scatter(dV_flex_order(:,4),TOF_flex_order(:,4),[],colors(12,:),'filled',...
    'DisplayName','Ast 4 Variations and Reordering'); % all changes on 4th ast
xlabel('$\Delta V$ [km/s]'); ylabel('TOF [d]'); 
legend('show','Location','northoutside','NumColumns',2)
