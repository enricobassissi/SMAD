%% --------------------------------------------------------------------- %%
%% ---------------------- Earth Ast1 Ast2 Transfer --------------------- %%
%% ------------------------- ARCH 1+4, LT RV --------------------------- %%
%% ------------------------------ SOO ---------------------------------- %%
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

sim.case_name = 'ARCH ID 6: LOW THRUST RENDEZVOUS ON EVERY ASTEROID + COASTING TIME';

%% add path of functions and python stuff
str_path=split(pwd, 'TrajOptimisation\LowThrust\Asteroids');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
str_path=split(pwd, 'Asteroids');
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
% load('data_elements_matrix_9.mat')
% data.p_number = 4;
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix, data.p_number);
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

load('data_processed_9_4.mat')

%% simulation parameters
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 200;                    % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 

sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.direction = 1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
sim.M = 100; % SC wet mass [kg]
sim.M_pods = 1; % mass of the payloads+landing stuff [kg]

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 500*3600*24/sim.TU; 
bound.TOF1_max = 2*365*3600*24/sim.TU; 
% TOF2 (3)
bound.TOF2_min = 0.5*365*3600*24/sim.TU; 
bound.TOF2_max = 2*365*3600*24/sim.TU; 
% N REV (4)
bound.N_REV1_min = 0; %0
bound.N_REV1_max = 1; %3
% N REV 2 (5)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 1; %3
% ID Permutation (6)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 
% bound.IDP_max = length(data.asteroid_names);
% C3 stuff
% Constraint on C3 Launcher (7)
sim.C3_max = 20; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% alpha, azimuth (8)
bound.alpha_min = -deg2rad(180);
bound.alpha_max = deg2rad(180);
% beta, elvation (9)
bound.beta_min = -deg2rad(180);
bound.beta_max = deg2rad(180);
% coasting time1 (10)
bound.CT1_min = 30*3600*24/sim.TU;
bound.CT1_max = 80*3600*24/sim.TU;
% coasting time2 (11)
bound.CT2_min = 30*3600*24/sim.TU;
bound.CT2_max = 80*3600*24/sim.TU;
% TOF3 (12)
bound.TOF3_min = 0.5*365*3600*24/sim.TU; %600
bound.TOF3_max = 2*365*3600*24/sim.TU; 
% N REV 3 (13)
bound.N_REV3_min = 0; %0
bound.N_REV3_max = 1; %3
% coasting time3 (14)
bound.CT3_min = 30*3600*24/sim.TU;
bound.CT3_max = 80*3600*24/sim.TU;
% TOF4 (15)
bound.TOF4_min = 0.2*365*3600*24/sim.TU; %600
bound.TOF4_max = 2*365*3600*24/sim.TU; 
% N REV 4 (16)
bound.N_REV4_min = 0; %0
bound.N_REV4_max = 1; %3

% x = [MJD0,TOF1,TOF2,NREV,NREV2,IDP,v_inf_magn,alpha,beta,CT1,CT2,TOF3,NREV3,
%      CT3,TOF4,NREV4]
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.N_REV1_min, bound.N_REV2_min, bound.IDP_min, ...
    bound.v_inf_magn_min, bound.alpha_min, bound.beta_min, ...
    bound.CT1_min, bound.CT2_min, bound.TOF3_min, bound.N_REV3_min, ...
    bound.CT3_min, bound.TOF4_min, bound.N_REV4_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.N_REV1_max, bound.N_REV2_max, bound.IDP_max, ...
    bound.v_inf_magn_max, bound.alpha_max, bound.beta_max, ...
    bound.CT1_max, bound.CT2_max, bound.TOF3_max, bound.N_REV3_max, ...
    bound.CT3_max, bound.TOF4_max, bound.N_REV4_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (4) NREV1, (5) NREV2, (6) IDP, (13) NREV3, (16) NREV4 
constr.IntCon = [4,5,6,13,16];

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
options.MaxGenerations = 50; % ideal 100

options.FunctionTolerance = 1e-6; %1e-9
options.MaxStallGenerations = ceil(options.MaxGenerations/10);

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
FitnessFunction = @(x) ff_1RL(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Building the solution structure
sol.asteroid_1 = data.PermutationMatrix(x(6),1);
sol.asteroid_2 = data.PermutationMatrix(x(6),2);
sol.asteroid_3 = data.PermutationMatrix(x(6),3);
sol.asteroid_4 = data.PermutationMatrix(x(6),4);

sol.departure_date = mjd20002date(x(1)*sim.TU/(3600*24));
sol.TOF1 = x(2)*sim.TU/86400;
sol.TOF2 = x(3)*sim.TU/86400;
sol.TOF3 = x(12)*sim.TU/86400;
sol.TOF4 = x(15)*sim.TU/86400;
sol.CT1 = x(10)*sim.TU/86400;
sol.CT2 = x(11)*sim.TU/86400;
sol.CT3 = x(14)*sim.TU/86400;
sol.mission_duration_d = sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3+sol.CT3+sol.TOF4;
sol.mission_duration_y = sol.mission_duration_d/365;
sol.arrival_date = mjd20002date(x(1)*sim.TU/(3600*24)+sol.mission_duration_d);

sol.mass_fraction = Fval;
sol.Nrev = [x(4), x(5), x(13), x(16)];
sol.v_inf_magn = x(7)*sim.DU/sim.TU;
sol.alpha = x(8);
sol.beta = x(9);

%% characteristic quantities plot
[output, r_encounter, v_encounter, sol] = plot_ff_ea_4ast_LT_soo_NLI(x,sim,data,sol);

figure()
subplot(5,1,1)
plot(output.t*sim.TU/86400,output.Thrust(:,1));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xline(output.t(sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(2*sim.n_sol)*sim.TU/86400,'LineWidth',2);
xline(output.t(3*sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(4*sim.n_sol)*sim.TU/86400,'LineWidth',2);
xline(output.t(5*sim.n_sol)*sim.TU/86400,'LineWidth',2)
xlabel('Time [days]')
ylabel('In-plane Thrust [N]')

subplot(5,1,2)
plot(output.t*sim.TU/86400,180/pi*output.Thrust(:,2));
xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t*sim.TU/86400,output.Thrust(:,3));
xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t*sim.TU/86400,sqrt(output.Thrust(:,1).^2 + output.Thrust(:,3).^2));
xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t*sim.TU/86400,output.m);
xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('Mass [kg]')

%% orbit plots
day1 = [2028 1 1 0 0 0];
day2 = [2032 1 1 0 0 0];

t1 = date2mjd2000(day1);
t2 = date2mjd2000(day2);
times = linspace(t1,t2,1000);

for i=1:length(times)
    % Orbit earth
    [kep1,ksun] = uplanet(times(i),3);
    [r1(i,1:3),v1(i,1:3)] = sv_from_coe(kep1,ksun);
    r1(i,1:3) = r1(i,1:3)/sim.DU;
    
    % Orbit mars
    [kep2,ksun] = uplanet(times(i),4);
    [r2(i,1:3),~] = sv_from_coe(kep2,ksun);
    r2(i,1:3) = r2(i,1:3)/sim.DU;
    
    % Orbit ast1
    [kep_ast_1] = uNEO2(times(i),sol.asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
    [r_ast1(i,1:3), v_ast1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
    r_ast1(i,1:3) = r_ast1(i,1:3)/sim.DU;
    
    % Orbit ast2
    [kep_ast_2] = uNEO2(times(i),sol.asteroid_2,data); % [km,-,rad,rad,rad,wrapped rad]
    [r_ast2(i,1:3), v_ast2] = sv_from_coe(kep_ast_2,ksun); % km, km/s
    r_ast2(i,1:3) = r_ast2(i,1:3)/sim.DU;
    
    % Orbit ast3
    [kep_ast_3] = uNEO2(times(i),sol.asteroid_3,data); % [km,-,rad,rad,rad,wrapped rad]
    [r_ast3(i,1:3), v_ast3] = sv_from_coe(kep_ast_3,ksun); % km, km/s
    r_ast3(i,1:3) = r_ast3(i,1:3)/sim.DU;
    
    % Orbit ast4
    [kep_ast_4] = uNEO2(times(i),sol.asteroid_4,data); % [km,-,rad,rad,rad,wrapped rad]
    [r_ast4(i,1:3), v_ast4] = sv_from_coe(kep_ast_4,ksun); % km, km/s
    r_ast4(i,1:3) = r_ast4(i,1:3)/sim.DU;
end

% transfer orbits
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

r_transf_orbit_3  = [output.r.leg3.*cos(output.theta.leg3), ...
    output.r.leg3.*sin(output.theta.leg3), output.z.leg3];
R_transf_orbit_3 = rotate_local2ecplitic(r_encounter.astD2,r_transf_orbit_3,sim.n_sol,output.Href.leg3);

r_transf_orbit_4  = [output.r.leg4.*cos(output.theta.leg4), ...
    output.r.leg4.*sin(output.theta.leg4), output.z.leg4];
R_transf_orbit_4 = rotate_local2ecplitic(r_encounter.astD3,r_transf_orbit_4,sim.n_sol,output.Href.leg4);

figure()
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj')
hold on
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt3 = plot3(R_transf_orbit_3(:,1),R_transf_orbit_3(:,2),R_transf_orbit_3(:,3),...
    'Color',colors(1,:));
hpt3.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt4 = plot3(R_transf_orbit_4(:,1),R_transf_orbit_4(:,2),R_transf_orbit_4(:,3),...
    'Color',colors(1,:));
hpt4.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astD2(1),r_encounter.astD2(2),r_encounter.astD2(3),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astA3(1),r_encounter.astA3(2),r_encounter.astA3(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astD3(1),r_encounter.astD3(2),r_encounter.astD3(3),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_3)
plot3(r_encounter.astA4(1),r_encounter.astA4(2),r_encounter.astA4(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_4)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% Earth
plot3(r1(:,1),r1(:,2),r1(:,3),...
    '--','Color',colors(8,:),'DisplayName','Earth'); % geocentric equatorial frame 
% Mars
plot3(r2(:,1),r2(:,2),r2(:,3),...
    '--','Color',colors(2,:),'DisplayName','Mars');
% Asteroid 1
plot3(r_ast1(:,1),r_ast1(:,2),r_ast1(:,3),...
    '--','Color',colors(3,:),'DisplayName',sol.asteroid_1);
% Asteroid 2
plot3(r_ast2(:,1),r_ast2(:,2),r_ast2(:,3),...
    '--','Color',colors(4,:),'DisplayName',sol.asteroid_2);
% Asteroid 3
plot3(r_ast3(:,1),r_ast3(:,2),r_ast3(:,3),...
    '--','Color',colors(5,:),'DisplayName',sol.asteroid_3);
% Asteroid 4
plot3(r_ast4(:,1),r_ast4(:,2),r_ast4(:,3),...
    '--','Color',colors(6,:),'DisplayName',sol.asteroid_4);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

