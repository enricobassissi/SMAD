%% --------------------------------------------------------------------- %%
%% ------------------ FLEXIBILITY ANALYSIS ----------------------------- %%
%% --------------------- OVER LAUNCH DATE ------------------------------ %%
%% ------------------ AND COMBINATION OF THE SEQUENCE ------------------ %%
%% ---------------------- BUT SAME ASTEROIDS --------------------------- %%
%% --------------------------------------------------------------------- %%
%% Description
% We have a solution optimised. Now we say, if the launch date is delayed,
% what happens to the mission? Which is the impact on the cost (dV and TOF)
% of the mission if we start some days after?
% The ECSS Margin Philosophy requires a flexibility of 3 weeks around the
% optimal launch date, half before the actual date and half after

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

%% add path of functions and python stuff
str_path=split(pwd, 'PostAnalysis\Flexibility');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Flexibility');
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
load('data_elements_matrix_44_63_2SC.mat')
load('160kg_dry_64mN.mat')

%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
% sim.PS.Isp = 4500/sim.TU;  % non-dimensional specific impulse % simone
sim.M1_end = 160; % SC wet mass [kg] %%
sim.M2_end = 160; % SC wet mass [kg] %%
sim.M_pods = 3.5; % mass of the payloads + landing stuff [kg] %%
sim.max_Available_Thrust = 0.02; % 5 [mN], BepiColombo is 250 mN but it's much bigger

%% Flexibility Parameters
% --- delay
% let's make arbitrary that the max delay time is the one of Rosetta
% mission, one of the most delayed mission 730 days
% regulation over launch widnow: 21 days
delay_time = 21;

% --- launch windows vector
% that now it becomes a fixed variable and not anymore a degree of freedom
% if i have problem for the launch, if i posticipate the launch is it a
% problem? meaning how does the costs (dV, TOF) changes?
N_div_points = 1; % how many points do you want? fine or coarse grid
lw_vector = linspace(sol.departure_mjd2000-delay_time/2,sol.departure_mjd2000+delay_time/2,delay_time/N_div_points);

% % --- Permutations among the 4 asteroid fixed
% % Number of possible combination of 4 asteroids among the ones in the list
% data.HowMany_flex = factorial(length(asteroid_sequence)) / factorial(length(asteroid_sequence) - 4);
% [data.PermutationMatrix_flex, ~] = permnUnique(asteroid_sequence, 4);

% asteroids selected
data.asteroid_1 = sol.asteroid_1;
data.asteroid_2 = sol.asteroid_2;
data.asteroid_a = sol.asteroid_a;
data.asteroid_b = sol.asteroid_b;

data.Nrev1 = sol.Nrev(1);
data.Nrev2 = sol.Nrev(2);
data.Nreva = sol.Nrev(3);
data.Nrevb = sol.Nrev(4);
 
%% Boundaries
% Departure dates
% it's now fixed by the flex analysis
lower_coeff = 0.8;
upper_coeff = 1.2;
% TOF1 (1)
bound_flex.TOF1_min = sol.TOF1*lower_coeff*86400/sim.TU; %+- 20% opt sol
bound_flex.TOF1_max = sol.TOF1*upper_coeff*86400/sim.TU; 
% TOF2 (2)
bound_flex.TOF2_min = sol.TOF2*lower_coeff*86400/sim.TU;
bound_flex.TOF2_max = sol.TOF2*upper_coeff*86400/sim.TU; 
% TOFa (3)
bound_flex.TOFa_min = sol.TOFa*lower_coeff*86400/sim.TU; %600
bound_flex.TOFa_max = sol.TOFa*upper_coeff*86400/sim.TU; 
% TOFb (4)
bound_flex.TOFb_min = sol.TOFb*lower_coeff*86400/sim.TU; %0.5*365
bound_flex.TOFb_max = sol.TOFb*upper_coeff*86400/sim.TU; 
% C3 stuff
% Constraint on C3 Launcher (5)
sim.C3_max = 30; % km^2/s^2
bound_flex.v_inf_magn_min = 0;
bound_flex.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (6)
bound_flex.az_min = -pi;
bound_flex.az_max = pi;
% elevation (7)
bound_flex.el_min = -pi/5;
bound_flex.el_max = pi/5;

% Coasting time 1 (8)
bound_flex.CT1_min = sol.CT1*lower_coeff*86400/sim.TU; % 30 days
bound_flex.CT1_max = sol.CT1*upper_coeff*86400/sim.TU; 
% Coasting time 2 (9)
bound_flex.CTa_min = sol.CTa*lower_coeff*86400/sim.TU; % 30 days
bound_flex.CTa_max = sol.CTa*upper_coeff*86400/sim.TU; 

% x = [TOF1,TOF2,TOF3,TOF4,v_inf_magn,az,el,CT1,CTa]
bound_flex.lb = [bound_flex.TOF1_min, bound_flex.TOF2_min, ...
    bound_flex.TOFa_min, bound_flex.TOFb_min, ...
    bound_flex.v_inf_magn_min, bound_flex.az_min, bound_flex.el_min, ...
    bound_flex.CT1_min, bound_flex.CTa_min]; % Lower bound
bound_flex.ub = [bound_flex.TOF1_max, bound_flex.TOF2_max, ...
    bound_flex.TOFa_max, bound_flex.TOFb_max, ...
    bound_flex.v_inf_magn_max, bound_flex.az_max, bound_flex.el_max, ...
    bound_flex.CT1_max, bound_flex.CTa_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints

%% Options
options = optimoptions(@gamultiobj);
% options.PlotFcn = @gaplotpareto;
options.Display = 'none';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
% options.HybridFcn = 'fgoalattain';

% ---- remember to change the integer position in the vector
% options.CreationFcn = @int_pop_2RL_moo_flex;
% options.MutationFcn = @int_mutation_2RL_moo_flex;
% options.CrossoverFcn = @int_crossoverarithmetic_2RL_moo_flex;

options.PopulationSize = 300; % ideal 1000
options.ParetoFraction = 0.6;
options.MaxGenerations = 30; % ideal 100

options.FunctionTolerance = 1e-9;
options.MaxStallGenerations = ceil(options.MaxGenerations/2);

% Parallel pool
% Open the parallel pool
par_pool = gcp; 
if isempty(par_pool)
    poolsize = 0;
else
    poolsize = par_pool.NumWorkers;
end

options.UseParallel = true;

numberOfVariables = length(bound_flex.ub); % Number of decision variables

%% MOO but looping over the dep date
% initialize
OBJ1 = zeros(length(lw_vector),1);
OBJ2 = zeros(length(lw_vector),1);
idx = 0;
tic
% loop over MJD0 for sensitivity analysis over departure date
for i = 1:length(lw_vector)
    
    idx = idx+1;
    fprintf('Cycle %d \n',idx);
    
    sim.MJD0_flexibility = lw_vector(i)*86400/sim.TU;
    FitnessFunction = @(x) ff_2RL_all_indietro_moo_LW_flex(x,sim,data); % Function handle to the fitness function
    [xx_flex,Fval_flex,exitFlag,Output,population,score] = gamultiobj(FitnessFunction,numberOfVariables,constr.A, ...
        constr.b,constr.Aeq,constr.beq,bound_flex.lb,bound_flex.ub,constr.nonlcon,options);
    FVAL_stored{i} = Fval_flex;
    XX_stored{i} = xx_flex;
    
    % find knee solution
    knee_sol_Fval = sqrt(Fval_flex(:,1).^2+(Fval_flex(:,2)./500).^2);
    idx_knee = find(min(knee_sol_Fval) == knee_sol_Fval);
    idx_knee = idx_knee(1);
    x_flex = xx_flex(idx_knee,:);
    X_stored{i} = x_flex;

    % save the result
    flex{i}.MJD0 = sim.MJD0_flexibility;
    flex{i}.dep_date = mjd20002date(sim.MJD0_flexibility)';
    flex{i}.TOF_tot_max_D = max(x_flex(1)+x_flex(2)+x_flex(8),x_flex(3)+x_flex(4)+x_flex(9));
    flex{i}.TOF_tot_max_Y = flex{i}.TOF_tot_max_D/365;
    flex{i}.end_of_mission_max_date = mjd20002date(sim.MJD0_flexibility+flex{i}.TOF_tot_max_D)';
    flex{i}.TOF1 = x_flex(1);
    flex{i}.TOF2 = x_flex(2);
    flex{i}.TOFa = x_flex(3);
    flex{i}.TOFb = x_flex(4);
    flex{i}.CT1 = x_flex(8);
    flex{i}.CTa = x_flex(9);
    flex{i}.v_inf_magn = x_flex(5)*sim.DU/sim.TU;
    flex{i}.az = x_flex(6);
    flex{i}.el = x_flex(7);
    flex{i}.az_deg = rad2deg(x_flex(6));
    flex{i}.el_deg = rad2deg(x_flex(7));
    flex{i}.obj_fun1 = Fval_flex(idx_knee,1);
    flex{i}.obj_fun2 = Fval_flex(idx_knee,2);
    
    plot_values = plot_ff_2RL_all_indietro_LW_flex(x_flex,sim,data);
    flex{i}.max_T = max([max(plot_values.T_1_magn),max(plot_values.T_2_magn),max(plot_values.T_a_magn),max(plot_values.T_b_magn)]);
    flex{i}.max_MF = max(plot_values.mass_fract_SC2,plot_values.mass_fract_SC2);
    
end
el_time=toc;
clearvars i idx

%% plot the result
N = length(lw_vector);
date_vect = zeros(N,1);
% N = 10;
for i=1:N
    
    OBJ1(i,1) = flex{i}.max_MF;
    OBJ1(i,1) = flex{i}.obj_fun1;
    if OBJ1(i,1) == 0 || OBJ1(i,1) > 0.4
        OBJ1(i,1) = NaN;
    end

    OBJ2(i,1) = flex{i}.max_T*1000; % mN
%     OBJ2(i,1) = flex{i}.obj_fun2; % mN
    if OBJ2(i,1) == 0 || OBJ2(i,1) > 800
        OBJ2(i,1) = NaN;
    end
    
    % plot date
    date_vect(i,1) = datenum(datetime(mjd20002date(lw_vector(i))));
    if date_vect(i,1) == datenum(datetime(mjd20002date(sol.departure_mjd2000)))
        OBJ1(i,1) = NaN;
        OBJ2(i,1) = NaN;
        date_vect(i,1) = NaN;
    end
end
clearvars i N

figure('Name','Flexibility over Departure Time')
yyaxis left
plot(date_vect,OBJ1,'*','Color',colors(1,:))
hold on
plot(datenum(datetime(mjd20002date(sol.departure_mjd2000))),max(sol.mass_fract_SC2, sol.mass_fract_SC1),'*','Color',colors(3,:))
ylabel('Mass Fraction')
set(gca,'ycolor',colors(1,:)) 
yyaxis right
plot(date_vect,OBJ2,'o','Color',colors(2,:))
plot(datenum(datetime(mjd20002date(sol.departure_mjd2000))),max(sol.max_T_SC1, sol.max_T_SC2)*1000,'o','Color',colors(3,:))
ylabel('Max Thrust [mN]')
set(gca,'ycolor',colors(2,:)) 
% xlabel('Departure MJD2000')
ax = gca;
ax.XTick=date_vect(1:4:end) ;
xtickangle(30)
datetick('x','yyyy mmm dd','keepticks')
xlim ([date_vect(1) date_vect(end)])
xline (datenum(datetime(mjd20002date(sol.departure_mjd2000))),'--','Color',colors(3,:),'LineWidth',2)