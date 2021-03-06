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

%% simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7           ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim)^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 200;                    % number of computational nodes
sim.x         = linspace(0,1,sim.n_sol)';   % 
sim.out_shape = 2;                  % out-of-plane shape, 2 = Conway-Wall shape
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW)

sim.vinf = 0;

sim.TOF_imposed_flag = 0;

%% Propulsive system parameters
sim.PS.Is = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 1000; % SC mass [kg]
sim.hp = 3; 
sim.kp = 3; %It is used just for sim.out_shape = 1;

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
sim.soo_lim.N_REV_max = 1; %3
% vinf_mag
sim.soo_lim.vinf_mag_min = 0.75;
sim.soo_lim.vinf_mag_max = 2;
% alpha
sim.soo_lim.alpha_min = -pi;
sim.soo_lim.alpha_max =  pi;
% beta
sim.soo_lim.beta_min = - pi;
sim.soo_lim.beta_max =  pi;


% % x = [MJD0,TOF,N_REV]
% sim.moo_bound.lb = [sim.moo_lim.mjd2000_ed, sim.moo_lim.TOF1_min, sim.moo_lim.N_REV_min]; % Lower bound
% sim.moo_bound.ub = [sim.moo_lim.mjd2000_ld, sim.moo_lim.TOF1_max, sim.moo_lim.N_REV_max]; % Upper bound
% x = [MJD0,TOF,N_REV,hp,vinf_mag,alpha,beta]
sim.soo_bound.lb = [sim.soo_lim.mjd2000_ed, sim.soo_lim.TOF1_min, sim.soo_lim.N_REV_min, sim.soo_lim.vinf_mag_min,sim.soo_lim.alpha_min,sim.soo_lim.beta_min]; % Lower bound
sim.soo_bound.ub = [sim.soo_lim.mjd2000_ld, sim.soo_lim.TOF1_max, sim.soo_lim.N_REV_max, sim.soo_lim.vinf_mag_max,sim.soo_lim.alpha_max,sim.soo_lim.beta_max]; % Upper bound

%% Constraints
sim.soo_constr.A = []; % linear inequality constraints
sim.soo_constr.b = []; % linear inequality constraints
sim.soo_constr.Aeq = []; % linear equality constraints
sim.soo_constr.beq = []; % linear equality constraints
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
options.MaxGenerations = 100; % ideal 100

options.FunctionTolerance = 1e-9; %1e-6
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

%% Build the soo
FitnessFunction = @(x) ff_ea_ma_LT_soo(x,sim); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,sim.soo_constr.A, ...
    sim.soo_constr.b,sim.soo_constr.Aeq,sim.soo_constr.beq,sim.soo_bound.lb,...
    sim.soo_bound.ub,sim.soo_constr.nonlcon,sim.soo_constr.IntCon,options);
el_time_min_pp = toc/60;


%% plot
[output, r1_true, r2_true] = plot_ff_ea_ma_LT_soo(x,sim);

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
%JD_departure = x(knee_sol,1);
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


r3  = [output.r.*cos(output.l) output.r.*sin(output.l) output.z];
R3 = rotate_local2ecplitic(r1_true,r3,sim.n_sol,output.h_ref_v);

figure()
plot3(R3(:,1),R3(:,2),R3(:,3))
axis equal
grid on
hold on
% Earth
hold on
plot3(r1(:,1),r1(:,2),r1(:,3),'--g'); % geocentric equatorial frame 
% Mars
hold on
plot3(r2(:,1),r2(:,2),r2(:,3),'--r');
plot3(0,0,0,'oy')

