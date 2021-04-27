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
%% simulation parameters
sim.mu    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU    = 149597870.7           ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 
sim.out_shape = 2;                  % out-of-plane shape, 2 = Conway-Wall shape
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW)

%% Propulsive system parameters
sim.PS.Is = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 1000; % SC mass [kg]
sim.hp = 3; 
sim.kp = 3;

%% Boundaries
% Departure dates
sim.moo_lim.date_ed = [2028, 1, 1, 0, 0, 0];
sim.moo_lim.date_ld =  [2031, 1, 1, 0, 0, 0];
sim.moo_lim.mjd2000_ed = date2mjd2000(sim.moo_lim.date_ed)*3600*24/sim.TU;
sim.moo_lim.mjd2000_ld = date2mjd2000(sim.moo_lim.date_ld)*3600*24/sim.TU;
% TOF1
sim.moo_lim.TOF1_min = 600*3600*24/sim.TU; 
sim.moo_lim.TOF1_max = 1000*3600*24/sim.TU; 
% N REV
sim.moo_lim.N_REV_min = -0.5;
sim.moo_lim.N_REV_max = 3.4999;

% x = [MJD0,TOF,N_REV]
sim.moo_bound.lb = [sim.moo_lim.mjd2000_ed, sim.moo_lim.TOF1_min, sim.moo_lim.N_REV_min]; % Lower bound
sim.moo_bound.ub = [sim.moo_lim.mjd2000_ld, sim.moo_lim.TOF1_max, sim.moo_lim.N_REV_max]; % Upper bound

%% Constraints
sim.moo_constr.A = []; % linear inequality constraints
sim.moo_constr.b = []; % linear inequality constraints
sim.moo_constr.Aeq = []; % linear equality constraints
sim.moo_constr.beq = []; % linear equality constraints

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

options.PopulationSize = 1000; % ideal 1000

options.ParetoFraction = 0.7;
options.MaxGenerations = 100; % ideal 100

options.FunctionTolerance = 1e-6;
options.MaxStallGenerations = 30;

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
FitnessFunction = @(x) ff_ea_ma_LT(x,sim); % Function handle to the fitness function
numberOfVariables = length(sim.moo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = gamultiobj(FitnessFunction,numberOfVariables,sim.moo_constr.A, ...
    sim.moo_constr.b,sim.moo_constr.Aeq,sim.moo_constr.beq,sim.moo_bound.lb,sim.moo_bound.ub,options);
el_time_min_pp = toc/60;

%% Find the knee solution
Fval(:,1) = Fval(:,1)*sim.TU/(3600*24);
[knee_idx, d] = find_knee_solution(Fval);

% Plot Pareto Plot
figure('Name','GA MO Pareto Plot')
title('Pareto Points in Parameter Space')
h_pp = plot(Fval(:,1),Fval(:,2),'o','Color',colors(1,:));
hold on
h_kpp = plot(Fval(knee_idx,1),Fval(knee_idx,2),'o','Color',colors(2,:));
xlabel('$Obj_1: \ TOF$ [d]')
ylabel('$Obj_2: \ mf frac$ [trattino in mezzo]')
legend([h_pp,h_kpp],'Sub-Optim Sol','Knee Sol')
clearvars h_pp h_kpp

%% plot
output = plot_ff_ea_ma_LT(x(knee_idx,:),sim);

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
day1 = [2021 1 1 0 0 0];
day2 = [2023 1 1 0 0 0];

t1 = date2mjd2000(day1);
t2 = date2mjd2000(day2);
times = linspace(t1,t2,1000);

for i=1:length(times)
    % Orbit 1
    [kep1,ksun] = uplanet(times(i),3);
    [r1(i,1:3),~] = sv_from_coe(kep1,ksun);
    r1(i,1:3) = r1(i,1:3)/sim.DU;
    % Orbit 2
    [kep2,ksun] = uplanet(times(i),4);
    [r2(i,1:3),~] = sv_from_coe(kep2,ksun);
    r2(i,1:3) = r2(i,1:3)/sim.DU;
end


figure()
plot3(output.r.*cos(output.l),output.r.*sin(output.l),output.z)
axis equal
grid on
hold on
% Earth
hold on
h(1)=plot3(r1(:,1),r1(:,2),r1(:,3),'--');
% Mars
hold on
h(2)=plot3(r2(:,1),r2(:,2),r2(:,3),'--');


% vedi in che ref frame è output.r 