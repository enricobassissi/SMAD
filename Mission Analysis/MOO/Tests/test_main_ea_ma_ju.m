%% ------------------------------------------ %%
%% ----------- EA MA JU Transfer ------------ %%
%% ------------------------------------------ %%
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
addpath time
addpath function
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
          30   51   120;... % TRUSTY AZURE +2
          0    103  98;... % PURE TEAL +2
          51   94   111;... % DEEP SPACE -1
          0    0    0]./255; % BLACK

%% INTRO
sim.mu = 1.32712440017987e11; % Sun planetary constant (mu = mass * G) (from DE405) [km^3/s^2]
sim.DU = 149597870.691; % Distance Unit = Astronomical Unit (AU) (from DE405) [km]
sim.TU = (sim.DU^3/sim.mu)^0.5; % Time Unit
sim.mu = 1;

%% Boundaries
% Departure dates
sim.moo_lim.date_ed = [2021, 1, 1, 0, 0, 0];
sim.moo_lim.date_ld =  [2025, 1, 1, 0, 0, 0];
sim.moo_lim.mjd2000_ed = date2mjd2000(sim.moo_lim.date_ed);
sim.moo_lim.mjd2000_ld = date2mjd2000(sim.moo_lim.date_ld);
% TOF1
sim.moo_lim.TOF1_min = 30; % days
sim.moo_lim.TOF1_max = 780; % days, 1/(1/365-1/686), one synodic period earth-mars
% Launcher velocity given and angles
sim.moo_lim.v_inf_magn_min = 0;
sim.moo_lim.v_inf_magn_max = sqrt(60); % c3 = 60 km/s^2
sim.moo_lim.alpha_min = deg2rad(0);
sim.moo_lim.alpha_max = deg2rad(360);
sim.moo_lim.beta_min = deg2rad(0);
sim.moo_lim.beta_max = deg2rad(360);
% Buffer time
sim.moo_lim.buffer_time_min = 60;
sim.moo_lim.buffer_time_max = 180;
% TOF2
sim.moo_lim.TOF2_min = 1*365; % days
sim.moo_lim.TOF2_max = 5*365; % days, 1/(1/365-1/4333), one synodic period earth-mars

% x = [MJD0, TOF1, v_inf_magn, alpha, beta, buffer_time, TOF2]
sim.moo_bound.lb = [sim.moo_lim.mjd2000_ed, sim.moo_lim.TOF1_min, sim.moo_lim.v_inf_magn_min,...
      sim.moo_lim.alpha_min, sim.moo_lim.beta_min, sim.moo_lim.buffer_time_min,...
      sim.moo_lim.TOF2_min]; % Lower bound
sim.moo_bound.ub = [sim.moo_lim.mjd2000_ld, sim.moo_lim.TOF1_max, sim.moo_lim.v_inf_magn_max,...
      sim.moo_lim.alpha_max, sim.moo_lim.beta_max, sim.moo_lim.buffer_time_max,...
      sim.moo_lim.TOF2_max]; % Upper bound

%% Constraints
sim.moo_constr.A=[]; % linear inequality constraints
sim.moo_constr.b= []; % linear inequality constraints
sim.moo_constr.Aeq = []; % linear equality constraints
sim.moo_constr.beq = []; % linear equality constraints

%% Options
options = optimoptions(@gamultiobj);
% options.PlotFcn = @gaplotpareto;
options.DistanceMeasureFcn = {@distancecrowding,'genotype'};
options.PopulationSize = 200;
options.ParetoFraction = 0.5;
options.MaxGenerations = 500;
options.FunctionTolerance = 1e-6;
options.MaxStallGenerations = 50;

%% Build the moo
FitnessFunction = @ff; % Function handle to the fitness function
numberOfVariables = length(sim.moo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = gamultiobj(FitnessFunction,numberOfVariables,sim.moo_constr.A, ...
    sim.moo_constr.b,sim.moo_constr.Aeq,sim.moo_constr.beq,sim.moo_bound.lb,sim.moo_bound.ub,options);
el_time_pp = toc;

% Find the knee solution
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
sol.MJD0 = x(knee_idx,1);
sol.dep_date = mjd20002date(sol.MJD0)';
sol.end_of_mission_date = mjd20002date(sol.MJD0+Fval(knee_idx,2))';
sol.dV_tot = Fval(knee_idx,1);
sol.TOF_tot = Fval(knee_idx,2);
sol.TOF1 = x(knee_idx,2);
sol.buffer_time = x(knee_idx,6);
sol.TOF2 = x(knee_idx,7);
sol.v_inf_magn = x(knee_idx,3);
sol.v_inf_alpha = rad2deg(x(knee_idx,4));
sol.v_inf_beta = rad2deg(x(knee_idx,5));

% retrieve info about that mission profile for table
knee_mission = [sol.dep_date(1:3);sol.end_of_mission_date(1:3);sol.dV_tot;sol.TOF_tot;sol.TOF1;...
                sol.buffer_time;sol.TOF2;sol.v_inf_magn;sol.v_inf_alpha;sol.v_inf_beta];
% print table values
Quantity = {'dep Y';'dep M';'dep D';'arr Y';'arr M';'arr D';'dV tot [km/s]';'TOF tot [d]';...
            'TOF1 [d]';'buffer time [d]';'TOF2 [d]';'|v inf| [km/s]';'v inf alpha [deg]';...
            'v inf beta [deg]'};
table(Quantity,knee_mission)

%% Plot trajectories
plot_mission(sol,colors)
