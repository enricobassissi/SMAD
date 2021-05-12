%% Initializing the Environment
clear; close all; clc;
% addpath time
% addpath function
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
sim.mu    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU    = 149597870.7           ; % distance unit [km]
% sim.mu    = muSun                   ; % actractor parameter [km^3 s^-2]
% sim.DU    = AU                      ; % distance unit [km]
sim.TU    = (sim.DU^3/sim.mu )^0.5; % time unit [s]
sim.mu    = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol = 100;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = -1;                     % direction of integration (1 FW, -1 BW)


str_path=split(pwd, 'TrajOptimisation\LowThrust\JP_LT\EARTH-NEREUS_SOO_NLI');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'EARTH-NEREUS_SOO_NLI');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
% try 
%     module = py.importlib.import_module('neo_api_function');
% catch
%     copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
%     module = py.importlib.import_module('neo_api_function');
% end

% [data_nereus,t_vector_nereus] = find_eph_neo("1982DB")
load('data_nereus.mat');

%% Propulsive system parameters
sim.PS.Isp = 3000/sim.TU;  % non-dimensional specific impulse

sim.M = 1000; % SC mass [kg]


%% Boundaries
% Departure dates
sim.soo_lim.date_ed = [2041, 1, 1, 0, 0, 0]; %30
sim.soo_lim.date_ld =  [2043, 1, 1, 0, 0, 0]; %50
sim.soo_lim.mjd2000_ed = date2mjd2000(sim.soo_lim.date_ed)*3600*24/sim.TU;
sim.soo_lim.mjd2000_ld = date2mjd2000(sim.soo_lim.date_ld)*3600*24/sim.TU;
% TOF1
sim.soo_lim.TOF1_min = 500*3600*24/sim.TU; %500
sim.soo_lim.TOF1_max = 900*3600*24/sim.TU; %1500
% N REV
sim.soo_lim.N_REV_min = 1; %0
sim.soo_lim.N_REV_max = 2;
% vinf_mag
sim.soo_lim.vinf_mag_min = 0*sim.TU/sim.DU;
sim.soo_lim.vinf_mag_max = 6*sim.TU/sim.DU;
% alpha
sim.soo_lim.alpha_min = -pi/2;
sim.soo_lim.alpha_max =  pi/2;
% beta
sim.soo_lim.beta_min = - pi/2;
sim.soo_lim.beta_max =  pi/2;



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
FitnessFunction = @(x) ff_soo_NEREUS(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(sim.soo_bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,sim.soo_constr.A, ...
    sim.soo_constr.b,sim.soo_constr.Aeq,sim.soo_constr.beq,sim.soo_bound.lb,...
    sim.soo_bound.ub,sim.soo_constr.nonlcon,sim.soo_constr.IntCon,options);
el_time_min_pp = toc/60;

dep_opt = mjd20002date(x(1)*sim.TU/(3600*24))
TOF_opt = x(2)*sim.TU/(3600*24) 
Nrev_opt = x(3)
vinf_opt = x(4)*sim.DU/sim.TU
alpha_opt = rad2deg(x(5))
beta_opt = rad2deg(x(6))
Fval 
%% plot
[output, r1_true, r2_true] = plot_ff_soo_NEREUS(x,sim,data);

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.Thrust(:,2));
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.Thrust(:,3));
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
xlabel('Time [days]')
ylabel('Mass [kg]')

massfrac_opt = Fval
TOF1_opt = x(2)*sim.TU/(3600*24)
%%
%JD_departure = x(knee_sol,1);
day1 = [2028 1 1 0 0 0];
day2 = [2030 1 1 0 0 0];

t1 = date2mjd2000(day1);
t2 = date2mjd2000(day2);
times = linspace(t1,t2,1000);

for i=1:length(times)
    % Orbit 1
    [kep1,ksun] = uplanet(times(i),3);
    [r1(i,1:3),v1(i,1:3)] = sv_from_coe(kep1,ksun);
    r1(i,1:3) = r1(i,1:3)/sim.DU;
    
    nereus_elements = interp1(data.t_vector, data.nereus_tbl, times(i), 'spline'); 
    [r2(i,1:3), v2] = sv_from_coe(nereus_elements,ksun); % km,km/s
    % adimensionalise
    r2(i,1:3) = r2(i,1:3)/sim.DU;
    
end


r3  = [output.r.*cos(output.theta) output.r.*sin(output.theta) output.z];
R3 = rotate_local2ecplitic(r1_true,r3,sim.n_sol,output.Href);

figure()
plot3(R3(:,1),R3(:,2),R3(:,3))
hold on
plot3(r1_true(1),r1_true(2),r1_true(3),'*m')
plot3(r2_true(1),r2_true(2),r2_true(3),'*c')
axis equal
grid on
hold on
% Earth
hold on
plot3(r1(:,1),r1(:,2),r1(:,3),'--g'); 
plot3(r2(:,1),r2(:,2),r2(:,3),'--b'); 

