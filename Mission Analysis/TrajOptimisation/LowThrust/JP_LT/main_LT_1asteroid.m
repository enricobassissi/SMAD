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
          30   51   120;... % (9) TRUSTY AZURE +2
          0    103  98;... % (10) PURE TEAL +2
          51   94   111;... % (11) DEEP SPACE -1
          0    0    0]./255; % (12) BLACK
      
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

str_path=split(pwd, 'TrajOptimisation\LowThrust\JP_LT');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));

%% simulation parameters
sim.mu    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW)

load('data.mat')
%% Propulsive system parameters
sim.PS.Is = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 100; % SC mass [kg]


%% Boundaries
% Departure dates(1)
sim.moo_lim.date_ed = [2022, 1, 1, 0, 0, 0];
sim.moo_lim.date_ld =  [2028, 1, 1, 0, 0, 0];
sim.moo_lim.mjd2000_ed = date2mjd2000(sim.moo_lim.date_ed)*3600*24/sim.TU;
sim.moo_lim.mjd2000_ld = date2mjd2000(sim.moo_lim.date_ld)*3600*24/sim.TU;
% TOF1(2)
sim.moo_lim.TOF1_min = 100*3600*24/sim.TU; 
sim.moo_lim.TOF1_max = 3*365*3600*24/sim.TU; 
% N REV(3)
sim.moo_lim.N_REV_min = 2;
sim.moo_lim.N_REV_max = 3;
% vinf_mag
sim.moo_lim.vinf_mag_min = 1.1/sim.DU*sim.TU;
sim.moo_lim.vinf_mag_max = sqrt(40)/sim.DU*sim.TU;
% alpha
sim.moo_lim.alpha_min = -pi;
sim.moo_lim.alpha_max =  pi;
% beta
sim.moo_lim.beta_min = - pi;
sim.moo_lim.beta_max =  pi;


% x = [MJD0,TOF,N_REV,vinf_mag,alpha,beta]
sim.moo_bound.lb = [sim.moo_lim.mjd2000_ed, sim.moo_lim.TOF1_min, sim.moo_lim.N_REV_min, sim.moo_lim.vinf_mag_min,sim.moo_lim.alpha_min,sim.moo_lim.beta_min]; % Lower bound
sim.moo_bound.ub = [sim.moo_lim.mjd2000_ld, sim.moo_lim.TOF1_max, sim.moo_lim.N_REV_max, sim.moo_lim.vinf_mag_max,sim.moo_lim.alpha_max,sim.moo_lim.beta_max]; % Upper bound

%% Constraints
sim.moo_constr.A = []; % linear inequality constraints
sim.moo_constr.b = []; % linear inequality constraints
sim.moo_constr.Aeq = []; % linear equality constraints
sim.moo_constr.beq = []; % linear equality constraints
sim.soo_constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
sim.soo_constr.IntCon = [3];

%% Options
options = optimoptions(@ga);
% options.PlotFcn = @gaplotpareto;
options.Display = 'iter';
% y = score; -> phenotype, Measure the distance in fitness function space; 
% y = pop; -> genotype, Measure the distance in decision variable space.
%options.DistanceMeasureFcn = {@distancecrowding,'phenotype'};
% A hybrid function is another minimization function that runs after the 
% multiobjective genetic algorithm terminates
% options.HybridFcn = 'fgoalattain';


options.PopulationSize = 1000; % ideal 1000
%options.ParetoFraction = 0.5;
options.MaxGenerations = 30; % ideal 100

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

%% Build the moo
FitnessFunction = @(x) ff_LT_1asteroid(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(sim.moo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,sim.moo_constr.A, ...
    sim.moo_constr.b,sim.moo_constr.Aeq,sim.moo_constr.beq,sim.moo_bound.lb,...
    sim.moo_bound.ub,sim.soo_constr.nonlcon,sim.soo_constr.IntCon,options);
el_time_min_pp = toc/60;


%% Find the knee solution
% Fval(:,1) = Fval(:,1)*sim.TU/(3600*24);
% [knee_idx, d] = find_knee_solution(Fval);


%% plot
[output, r1_true, r2_true,v1_true,v2_true] = plot_LT_1asteroid(x,sim,data);
output.u = output.Thrust;
output.l = output.theta;

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.u(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.u(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.u(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.u(:,1).^2 + output.u(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')


%%

r3 = [output.r.*cos(output.l) output.r.*sin(output.l) output.z]; 


figure()
plot3(r3(:,1),r3(:,2),r3(:,3))
axis equal
grid on
hold on
plot3(r1_true(1),r1_true(2),r1_true(3),'m*')
plot3(r2_true(1),r2_true(2),r2_true(3),'c*')
