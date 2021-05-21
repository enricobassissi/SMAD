%% --------------------------------------------------------------------- %%
%% ------------------------- Earth Mars Transfer ----------------------- %%
%% ----------------------------- MOO ----------------------------------- %%
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
str_path=split(pwd, 'TrajOptimisation\LowThrust\EARTH-MARS_MOO_NLI');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'EARTH-MARS_MOO_NLI');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu_dim)^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW)

sim.vinf = 0; % Parabolic escape

sim.TOF_imposed_flag = 1;

%% Propulsive system parameters
sim.PS.Isp = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 1000; % SC mass [kg]

%% Boundaries
% Departure dates
sim.soo_lim.date_ed = [2028, 1, 1, 0, 0, 0];
sim.soo_lim.date_ld =  [2031, 1, 1, 0, 0, 0];
sim.soo_lim.mjd2000_ed = date2mjd2000(sim.soo_lim.date_ed)*3600*24/sim.TU;
sim.soo_lim.mjd2000_ld = date2mjd2000(sim.soo_lim.date_ld)*3600*24/sim.TU;
% TOF1
sim.soo_lim.TOF1_min = 600*3600*24/sim.TU; 
sim.soo_lim.TOF1_max = 1000*3600*24/sim.TU; 
% N REV
sim.soo_lim.N_REV_min = 0;
sim.soo_lim.N_REV_max = 2;

% % x = [MJD0,TOF,N_REV]
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed, sim.soo_lim.TOF1_min, sim.soo_lim.N_REV_min]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max, sim.soo_lim.N_REV_max]; % Upper bound

%% Constraints
sim.soo_constr.A = []; % linear inequality constraints
sim.soo_constr.b = []; % linear inequality constraints
sim.soo_constr.Aeq = []; % linear equality constraints
sim.soo_constr.beq = []; % linear equality constraints
sim.soo_constr.nonlcon = []; % linear equality constraints

%% Options
options = optimoptions(@gamultiobj);
% options.PlotFcn = @gaplotpareto;
options.Display = 'iter';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
options.HybridFcn = 'fgoalattain';
options.CreationFcn = @int_pop;
options.MutationFcn = @int_mutation;
options.CrossoverFcn = @int_crossoverarithmetic;

options.PopulationSize = 600; % ideal 1000
options.ParetoFraction = 0.5;
options.MaxGenerations = 30; % ideal 100

options.FunctionTolerance = 1e-6; %1e-9
options.MaxStallGenerations = 10;

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
FitnessFunction = @(x) ff_ea_ma_moo(x,sim); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output,population,score] = gamultiobj(FitnessFunction,numberOfVariables,sim.soo_constr.A, ...
    sim.soo_constr.b,sim.soo_constr.Aeq,sim.soo_constr.beq,sim.soo_bound.lb,...
    sim.soo_bound.ub,sim.soo_constr.nonlcon,options);
el_time_min_pp = toc/60;

%% Find pareto plot and knee solution
[knee_idx1] = find_knee_solution(Fval);
knee_idx = knee_idx1(1);
% knee_idx = find(min(Fval(:,1))==Fval(:,1));
adim2days = sim.TU/(3600*24);

% Plot Pareto Plot
figure('Name','GA MO Pareto Plot')
title('Pareto Points in Parameter Space')
h_pp = plot(Fval(:,1),Fval(:,2)*adim2days,'o','Color',colors(1,:));
hold on
h_kpp = plot(Fval(knee_idx,1),Fval(knee_idx,2)*adim2days,'o','Color',colors(2,:));
xlabel('$Obj_1: \ mass_{frac}$ [km/s]')
ylabel('$Obj_2: \ TOF$ [d]')
legend([h_pp,h_kpp],'Sub-Optim Sol','Knee Sol')
clearvars h_pp h_kpp

% display final objetive values
massfrac_opt = Fval(knee_idx,1)
TOF1_opt = Fval(knee_idx,2)*adim2days
N_REV_chosen = x(knee_idx,3)
%% plot
[output, r1_true, r2_true] = plot_ff_ea_ma_moo(x(knee_idx,:),sim);

figure()
subplot(3,2,1)
plot(output.t*adim2days,output.Thrust(:,1));
xlabel('Time [days]')
ylabel('In-plane T [N]')

subplot(3,2,2)
plot(output.t*adim2days,rad2deg(output.Thrust(:,2)));
xlabel('Time [days]')
ylabel('In-plane T angle [deg]')

subplot(3,2,3)
plot(output.t*adim2days,output.Thrust(:,3));
xlabel('Time [days]')
ylabel('out-of-plane T [N]')

subplot(3,2,4)
plot(output.t*adim2days,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
xlabel('Time [days]')
ylabel('T [N]')

subplot(3,2,5)
plot(output.t*adim2days,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

%% ORBIT PLOT
% complete planet orbits
day1 = [2028 1 1 0 0 0];
day2 = [2031 1 1 0 0 0];
t1 = date2mjd2000(day1);
t2 = date2mjd2000(day2);
times = linspace(t1,t2,1000);
for i=1:length(times)
    % Orbit 1
    [kep1,ksun] = uplanet(times(i),3);
    [r1(i,1:3),v1(i,1:3)] = sv_from_coe(kep1,ksun);
    r1(i,1:3) = r1(i,1:3)/sim.DU;
    % Orbit 2
    [kep2,ksun] = uplanet(times(i),4);
    [r2(i,1:3),~] = sv_from_coe(kep2,ksun);
    r2(i,1:3) = r2(i,1:3)/sim.DU;
end

% Position at moment of encounter
mjd2000_dep = x(knee_idx,1)*adim2days;
date_dep = mjd20002date(mjd2000_dep);
mjd2000_arr = x(knee_idx,1)*adim2days + TOF1_opt;
date_arr = mjd20002date(mjd2000_arr);
% earth
[kep_dep_ea,ksun] = uplanet(mjd2000_dep,3);
[r_dep_ea(1:3),v_dep_ea(1:3)] = sv_from_coe(kep_dep_ea,ksun);
r_dep_ea(1:3) = r_dep_ea(1:3)/sim.DU;
% mars
[kep_arr_ma,ksun] = uplanet(mjd2000_arr,4);
[r_arr_ma(1:3),v_arr_ma(1:3)] = sv_from_coe(kep_arr_ma,ksun);
r_arr_ma(1:3) = r_arr_ma(1:3)/sim.DU;

% trajectory
r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
R3 = rotate_local2ecplitic(r1_true,r3,sim.n_sol,output.Href);

% plots
figure()
plot3(R3(:,1),R3(:,2),R3(:,3),'Color',colors(1,:),'DisplayName','Trajectory')
axis equal
grid on
hold on
% Earth
plot3(r1(:,1),r1(:,2),r1(:,3),'--','Color',colors(3,:),'DisplayName','Earth'); % heliocentric ecliptic frame 
plot3(r_dep_ea(:,1),r_dep_ea(:,2),r_dep_ea(:,3),'^','Color',colors(3,:),'DisplayName','Dep Point');
hp2 = plot3(r1_true(1),r1_true(2),r1_true(3),'o','Color',colors(3,:));
hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Mars
plot3(r2(:,1),r2(:,2),r2(:,3),'--','Color',colors(2,:),'DisplayName','Mars');
plot3(r_arr_ma(:,1),r_arr_ma(:,2),r_arr_ma(:,3),'^','Color',colors(2,:),'DisplayName','Arr Point');
hp3 = plot3(r2_true(1),r2_true(2),r2_true(3),'o','Color',colors(2,:));
hp3.Annotation.LegendInformation.IconDisplayStyle = 'off';
% sun
hsun = plot3(0,0,0,'o','Color',colors(4,:));
hsun.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('show')
