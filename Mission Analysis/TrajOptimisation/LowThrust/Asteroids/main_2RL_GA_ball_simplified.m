%% --------------------------------------------------------------------- %%
%% ------------------------ Earth Ast1 Ast2  --------------------------- %%
%% ------------------------ Earth Asta Astb ---------------------------- %%
%% ------------------------- ARCH 2SC, LT RV --------------------------- %%
%% ------------------------------ SOO ---------------------------------- %%
%% - UNPOWERED THRUST LEG FROM EARTH DEPARTURE AND GA AGAIN WITH EARTH - %%
%% ---- only one sc analysed ---- %%
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

sim.case_name = 'ARCH ID 8: LOW THRUST FLYBY 2 SC, EACH ON TWO ASTEROIDs';

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
% load('data_elements_matrix_42_61_2SC.mat')
% data.p_number = 2;
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix, data.p_number);
% 
% data.data_element_matrix = data_elements_matrix;
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

load('data_processed_42_61_2SC.mat')
load('power_propulsion_data.mat')
power_propulsion_data.T_limit = 20; % mN

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
sim.ID_FLYBY = 3;
sim.TOF_imposed_flag = 1;
sim.PS.Isp = 3200/sim.TU;  % non-dimensional specific impulse
sim.M1 = 100; % SC wet mass [kg] %%
% sim.M2 = 100; % SC wet mass [kg] %%
sim.M_pods = 5; % mass of the payloads + landing stuff [kg] + propellant
sim.max_Available_Thrust = 0.015; % 10 [mN], BepiColombo is 250 mN but it's much bigger
sim.dV_launch_man = 0.08/sim.DU*sim.TU; % 80 m/s max manouvre at beginning 
sim.C3_max = 30/(sim.DU^2)*(sim.TU^2); % max C3 launcher, km^2/s^2, adimensionalised
sim.dV_DSM_max = 0.1*sim.TU/sim.DU; % max 200 m/s

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;

% -- SC1 -- %
% DSM stuff SC1, ballistic
% TOF_DSM1 (2), from earth dep to DSM time
bound.TOF_DSM1_min = 0; % days
bound.TOF_DSM1_max = 1*365*3600*24/sim.TU; % days
% rD (3)
bound.rD1_min = 0.8; % AU, +/- 20%
bound.rD1_max = 1.2; % AU
% thetaD (4)
bound.thetaD1_min = 0; % rad
bound.thetaD1_max = 2*pi; % rad

% GA Stuff SC1, ballistic
% TOFGA (5), from DSM1 to GA 
bound.TOFGA_1_min = 0.1*365*3600*24/sim.TU; %10*3600
bound.TOFGA_1_max = 1.2*365*3600*24/sim.TU; 
% azimuth out (POST-GRAVITY ASSIST)(6)
bound.az_out_GA_1_min = -pi;
bound.az_out_GA_1_max =  pi;
% elevation (POST-GRAVITY ASSIST)(7)
bound.el_out_GA_1_min = -pi/2;
bound.el_out_GA_1_max =  pi/2;

% From gravity assist SC1 to ast 1, start of LT powering
% TOF1 (8)
bound.TOF1_min = 0.2*365*3600*24/sim.TU; %0.2 
bound.TOF1_max = 3*365*3600*24/sim.TU; 
% N REV 1 (9)
bound.N_REV1_min = 0; %0
bound.N_REV1_max = 2; %3
% Coasting time 1 (10)
bound.CT1_min = 30*3600*24/sim.TU; % 30 days
bound.CT1_max = 80*3600*24/sim.TU; 
% From ast1+CT1 to ast 2
% TOF2 (11)
bound.TOF2_min = 0.2*365*3600*24/sim.TU; %0.2
bound.TOF2_max = 3*365*3600*24/sim.TU; 
% N REV 2 (12)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 2; %3

% -- ID Permutation (13)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 

% x = [
% Departure dates (1)
% TOF_DSM1 (2), from earth dep to DSM time
% rD (3)
% thetaD (4)
% TOFGA (5), from DSM1 to GA
% azimuth out (POST-GRAVITY ASSIST)(6)
% elevation (POST-GRAVITY ASSIST)(7)
% TOF1 (8)
% N REV 1 (9)
% Coasting time 1 (10)
% TOF2 (11)
% N REV 2 (12)

% -- ID Permutation (13)]

% for now we fix asteroids!!!
% bound.lb = [bound.mjd2000_ed,bound.TOF_DSM1_min,bound.rD1_min,bound.thetaD1_min,bound.TOFGA_1_min,...
%             bound.az_out_GA_1_min,bound.el_out_GA_1_min,bound.TOF1_min,...
%             bound.N_REV1_min,bound.CT1_min,bound.TOF2_min,bound.N_REV2_min,bound.IDP_min];
% bound.ub = [bound.mjd2000_ld,bound.TOF_DSM1_max,bound.rD1_max,bound.thetaD1_max,bound.TOFGA_1_max,...
%             bound.az_out_GA_1_max,bound.el_out_GA_1_max,bound.TOF1_max,...
%             bound.N_REV1_max,bound.CT1_max,bound.TOF2_max,bound.N_REV2_max,bound.IDP_max];
bound.lb = [bound.mjd2000_ed,bound.TOF_DSM1_min,bound.rD1_min,bound.thetaD1_min,bound.TOFGA_1_min,...
            bound.az_out_GA_1_min,bound.el_out_GA_1_min,bound.TOF1_min,...
            bound.N_REV1_min,bound.CT1_min,bound.TOF2_min,bound.N_REV2_min];
bound.ub = [bound.mjd2000_ld,bound.TOF_DSM1_max,bound.rD1_max,bound.thetaD1_max,bound.TOFGA_1_max,...
            bound.az_out_GA_1_max,bound.el_out_GA_1_max,bound.TOF1_max,...
            bound.N_REV1_max,bound.CT1_max,bound.TOF2_max,bound.N_REV2_max];

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (9) NREV1, (12) NREV2, (13) IDP
% constr.IntCon = [9,12,13];
constr.IntCon = [9,12];

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
options.MaxGenerations = 300; % ideal 100

options.FunctionTolerance = 1e-6; %1e-9
options.MaxStallGenerations = ceil(options.MaxGenerations/15);

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
FitnessFunction = @(x) ff_2RL_GA_ball_simplified(x,sim,data,power_propulsion_data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Building the solution structure
sol.departure_date = mjd20002date(x(1)*sim.TU/86400);
sol.departure_mjd2000 = x(1)*sim.TU/86400;
% 1st spacecraft characteristic times
sol.TOF_DSM1 = x(2)*sim.TU/86400; % It's the time from Earth to DSM for SC1
sol.TOFGA_1 = x(5)*sim.TU/86400; % DSM to GA
sol.TOF1 = x(8)*sim.TU/86400; % tof sc1 GA to Ast1 
sol.CT1 = x(10)*sim.TU/86400;
sol.TOF2 = x(11)*sim.TU/86400; % tof sc1 from ast 1 to 2nd asteroid
% mission duration
sol.mission_duration_d_SC1 = sol.TOF_DSM1+sol.TOFGA_1+sol.TOF1+sol.CT1+sol.TOF2;
sol.mission_duration_y_SC1 = sol.mission_duration_d_SC1/365;
sol.arrival_dat = mjd20002date(x(1)*sim.TU/(3600*24)+sol.mission_duration_d_SC1);
% -- DSM SC1
sol.rD_mag1 = x(3); % AU
sol.thetaD1 = x(4); % rad
% asteorids
% sol.asteroid_1 = data.PermutationMatrix(x(13),1);
% sol.asteroid_2 = data.PermutationMatrix(x(13),2);
sol.asteroid_1 = "2020VV";
sol.asteroid_2 = "2010UH";

sol.Nrev = [x(9), x(12)];

% [sol_dates] = sol_to_dates_of_mission_LT(sol,'2RL_GA_DSM')

%% characteristic quantities plot
[output, r_encounter, v_encounter, sol] = plot_ff_2RL_GA_ball_simplified(x,sim,data,sol,power_propulsion_data);

%% Thrust
sol.max_T_SC1 = max(output.T_magn_SC1);

figure()
subplot(5,1,1)
hp1 = plot(output.t_SC1*sim.TU/86400,output.Thrust_SC1(:,1),'Color',colors(1,:));
% hold on
% hp2 = plot(output.t_SC2*sim.TU/86400,output.Thrust_SC2(:,1),'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xline(output.t(sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(2*sim.n_sol)*sim.TU/86400,'LineWidth',2);
% xline(output.t(3*sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(4*sim.n_sol)*sim.TU/86400,'LineWidth',2);
% xline(output.t(5*sim.n_sol)*sim.TU/86400,'LineWidth',2)
% xlabel('Time [days]')
ylabel('In-plane Thrust [N]')
lgd = legend(hp1,'SC1');
lgd.Location = 'eastoutside';
lgd.NumColumns = 1;
clearvars lgd hp1 hp2

subplot(5,1,2)
plot(output.t_SC1*sim.TU/86400,output.Thrust_SC1(:,2),'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.Thrust_SC2(:,2),'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t_SC1*sim.TU/86400,output.Thrust_SC1(:,3),'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.Thrust_SC2(:,3),'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t_SC1*sim.TU/86400,output.T_magn_SC1,'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.T_magn_SC2,'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t_SC1*sim.TU/86400,output.m_SC1,'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.m_SC2,'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('Mass [kg]')

%% only mass e thrust
subplot(2,1,1)
plot(output.t_SC1*sim.TU/86400,output.m_SC1,'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.m_SC2,'Color',colors(2,:));
xline(sol.TOF1,'LineWidth',2,'LineStyle',':','Color',colors(1,:)); 
xline(sol.TOF1+sol.CT1,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
% xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(2,:)); 
% xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
% xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
xlabel('Time [days]')
ylabel('Mass [kg]')

subplot(2,1,2)
hp1 = plot(output.t_SC1*sim.TU/86400,output.T_magn_SC1,'Color',colors(1,:));
xline(sol.TOF1,'LineWidth',2,'LineStyle',':','Color',colors(1,:)); 
xline(sol.TOF1+sol.CT1,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
% xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(2,:)); 
% xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
% xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
xlabel('Time [days]')
ylabel('Thrust [N]')

lgd = legend(hp1,'SC1');
% lgd.Location = 'eastoutside';
lgd.NumColumns = 1;
clearvars lgd hp1

%% orbit plots
% transfer orbits
% SC1
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.GA1,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

% SC 1-First leg: Earth -> DSM
departure_earth_sec = sol.departure_mjd2000*3600*24;
arrival_dsm_sec = (sol.departure_mjd2000+sol.TOF_DSM1)*3600*24;
y0_EA_DSM1 = [r_encounter.EA*sim.DU; v_encounter.dep_EA_DSM1*sim.DU/sim.TU]; %km, km/s; velocity from lambert arc transfer orbit injection
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_EA_DSM,y_EA_DSM] = ode113(@rates, [departure_earth_sec arrival_dsm_sec], y0_EA_DSM1,options,'sun');
% SC 1-First leg: Earth -> DSM
departure_dsm_sec = (sol.departure_mjd2000+sol.TOF_DSM1)*3600*24;
arrival_GA_sec = (sol.departure_mjd2000+sol.TOF_DSM1+sol.TOFGA_1)*3600*24;
y0_DSM1_GA = [r_encounter.DSM1*sim.DU; v_encounter.dep_DSM1_GA1*sim.DU/sim.TU]; %km, km/s; velocity from lambert arc transfer orbit injection
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_DSM_GA,y_DSM_GA] = ode113(@rates, [departure_dsm_sec arrival_GA_sec], y0_DSM1_GA,options,'sun');

figure()
% plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
%     'Color',colors(1,:),'DisplayName','Traj SC1')
hold on
% plot3( y_EA_DSM(:,1)./sim.DU, y_EA_DSM(:,2)./sim.DU, y_EA_DSM(:,3)./sim.DU,'Color',colors(2,:),...
%     'DisplayName','EA to DSM');
% hdsm1 = plot3(y_DSM_GA(:,1)./sim.DU, y_DSM_GA(:,2)./sim.DU, y_DSM_GA(:,3)./sim.DU,...
%     'Color',colors(3,:),'DisplayName','DSM to GA');
% hc1 = plot3( yC1(:,1)/sim.DU, yC1(:,2)/sim.DU, yC1(:,3)/sim.DU,'*','Color',colors(1,:),...
%     'Markersize',3,'DisplayName','CT');
% hc1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    'o','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.DSM1(1),r_encounter.DSM1(2),r_encounter.DSM1(3),...
    '^','Color',colors(5,:),'DisplayName','DSM')
plot3(r_encounter.GA1(1),r_encounter.GA1(2),r_encounter.GA1(3),...
    'd','Color',colors(8,:),'DisplayName','Earth GA')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1+' Arr')
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1+' Dep')
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2+' Arr')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+output.t1(end)*sim.TU/(3600*24);
hello_orbit2 = sol.departure_mjd2000+(output.t1(end)+ sol.CT1 + output.t2(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

%% post analysis
Traj_SC1 = [R_transf_orbit_1;R_transf_orbit_2];
max_dist_sun_SC1 = max(vecnorm(Traj_SC1,2,2));
min_dist_sun_SC1 = min(vecnorm(Traj_SC1,2,2));

Traj_SC2 = [R_transf_orbit_a;R_transf_orbit_b];
max_dist_sun_SC2 = max(vecnorm(Traj_SC2,2,2));
min_dist_sun_SC2 = min(vecnorm(Traj_SC2,2,2));

mjd_earth_plot_SC1 = [output.t1;output.t2].*sim.TU/86400 + sol.departure_mjd2000;
mjd_earth_plot_SC2 = [output.ta;output.tb].*sim.TU/86400 + sol.departure_mjd2000;
for k=1:length(mjd_earth_plot_SC1)
    [kep_SC1,ksun] = uplanet(mjd_earth_plot_SC1(k), 3);
    [r__E_SC1, ~] = sv_from_coe(kep_SC1,ksun);  
    R__E_SC1(k,:)=r__E_SC1/sim.DU; % it's in km it becomes AU, to be plotted
    
    [kep_SC2,ksun] = uplanet(mjd_earth_plot_SC2(k), 3);
    [r__E_SC2, ~] = sv_from_coe(kep_SC2,ksun);  
    R__E_SC2(k,:)=r__E_SC2/sim.DU; % it's in km it becomes AU, to be plotted
end
clearvars k r__E_SC1 r__E_SC2

dist_sc_earth_SC1 = Traj_SC1 - R__E_SC1;
max_dist_earth_SC1 = max(vecnorm(dist_sc_earth_SC1,2,2));
min_dist_earth_SC1 = min(vecnorm(dist_sc_earth_SC1,2,2));

dist_sc_earth_SC2 = Traj_SC2 - R__E_SC2;
max_dist_earth_SC2 = max(vecnorm(dist_sc_earth_SC2,2,2));
min_dist_earth_SC2 = min(vecnorm(dist_sc_earth_SC2,2,2));

sol.TOF1
output.t1(end)
transf1_check = abs(sol.TOF1- output.t1(end)*sim.TU/(3600*24))

sol.TOF2
output.t2(end)
transf2_check = abs(sol.TOF2- output.t2(end)*sim.TU/(3600*24))

sol.TOFa
output.ta(end)
transfa_check = abs(sol.TOFa - output.ta(end)*sim.TU/(3600*24))

sol.TOFb
output.tb(end)
transfb_check = abs(sol.TOFb - output.tb(end)*sim.TU/(3600*24))

% Propulsion for NIKITA
Propulsion.Magnitude_Thrust = [output.T_magn_SC1 output.T_magn_SC2];
Propulsion.InPlane_Thrust   = [output.Thrust_SC1(:,1) output.Thrust_SC2(:,1)];
Propulsion.gamma_angle      = [output.Thrust_SC1(:,2) output.Thrust_SC2(:,2)];
Propulsion.OutofPlane_Thrust = [output.Thrust_SC1(:,3) output.Thrust_SC2(:,3)];

    %% Extract the Sun-SpaceCraft Distance for all the trajectory [km]
%     sol.SunSpacecraftDistanceNorm = vecnorm([yEA1;y_ast12;y_EAa;y_astab],2,2); % 2,2 means norm 2 and by row
    sol.Spacecraft1Trajectory = [y_EA_DSM;yC1;y_ast12];
    sol.Spacecraft2Trajectory = [y_EAa;yCa;y_astab];
    sol.SC1time = [t_EA_DSM;tC1;t_ast12;];
    sol.SC2time = [t_EAa;tCa;t_astab];
%     [sol.angles.SAA,sol.angles.EVA,sol.angles.SCA,sol.angles.SolarConjunction] = aspect_angles(sol);

%% Plot with thrust vectors
Tlocal_transf_orbit_1  = [-sol.T_1(:,1).*sin(output.theta.leg1), ...
    sol.T_1(:,1).*cos(output.theta.leg1), sol.T_1(:,3)];
T_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit_1,sim.n_sol,output.Href.leg1);

Tlocal_transf_orbit_2  = [-sol.T_2(:,1).*sin(output.theta.leg2), ...
    sol.T_2(:,1).*cos(output.theta.leg2), sol.T_2(:,3)];
T_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,Tlocal_transf_orbit_2,sim.n_sol,output.Href.leg2);

Tlocal_transf_orbit_a  = [-sol.T_a(:,1).*sin(output.theta.lega), ...
    sol.T_a(:,1).*cos(output.theta.lega), sol.T_a(:,3)];
T_transf_orbit_a = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit_a,sim.n_sol,output.Href.lega);

Tlocal_transf_orbit_b  = [-sol.T_b(:,1).*sin(output.theta.legb), ...
    sol.T_b(:,1).*cos(output.theta.legb), sol.T_b(:,3)];
T_transf_orbit_b = rotate_local2ecplitic(r_encounter.astDa,Tlocal_transf_orbit_b,sim.n_sol,output.Href.legb);

figure('Name','Thrust Plot')
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1')
hold on
hpt1 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),...
    'Color',colors(2,:),'DisplayName','Traj SC2')
hpt2 = plot3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),...
    'Color',colors(2,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.astAa(1),r_encounter.astAa(2),r_encounter.astAa(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_a)
plot3(r_encounter.astDa(1),r_encounter.astDa(2),r_encounter.astDa(3),...
    '*','Color',colors(4,:),'DisplayName',sol.asteroid_a)
plot3(r_encounter.astAb(1),r_encounter.astAb(2),r_encounter.astAb(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_b)

quiver3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3),2)
quiver3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3),2)
quiver3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),T_transf_orbit_a(:,1),T_transf_orbit_a(:,2),T_transf_orbit_a(:,3),2)
quiver3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),T_transf_orbit_b(:,1),T_transf_orbit_b(:,2),T_transf_orbit_b(:,3),2)

axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
title('In-plane + out-of-plane Thrust')
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
%view(2)

figure('Name','Plane thrust 1')
hold on
plot(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),...
    'Color',colors(1,:),'DisplayName','Traj')
hpt1 = plot(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),...
    'Color',colors(1,:));
hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';

% plot(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),...
%     'Color',colors(1,:),'DisplayName','Traj')
% hpt2 = plot(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),...
%     'Color',colors(1,:));
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(r_encounter.EA(1),r_encounter.EA(2),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot(r_encounter.astA1(1),r_encounter.astA1(2),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot(r_encounter.astD1(1),r_encounter.astD1(2),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot(r_encounter.astA2(1),r_encounter.astA2(2),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
% plot(r_encounter.astAa(1),r_encounter.astAa(2),...
%     '^','Color',colors(5,:),'DisplayName',sol.asteroid_a)
% plot(r_encounter.astDa(1),r_encounter.astDa(2),...
%     '*','Color',colors(5,:),'DisplayName',sol.asteroid_a)
% plot(r_encounter.astAb(1),r_encounter.astAb(2),...
%     '^','Color',colors(6,:),'DisplayName',sol.asteroid_b)

quiver(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),2,...
    'Color',colors(3,:),'DisplayName',sol.asteroid_1)
quiver(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),2,...
    'Color',colors(4,:),'DisplayName',sol.asteroid_2)
% quiver(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),T_transf_orbit_a(:,1),T_transf_orbit_a(:,2),2,...
%     'Color',colors(5,:),'DisplayName',sol.asteroid_a)
% quiver(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),T_transf_orbit_b(:,1),T_transf_orbit_b(:,2),2,...
%     'Color',colors(6,:),'DisplayName',sol.asteroid_b)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]');
% title('In-plane Thrust(shall be tangential)')
% Sun
plot(0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

figure('Name','Plane thrust 2')
hold on
% plot(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),...
%     'Color',colors(1,:),'DisplayName','Traj')
% hpt1 = plot(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),...
%     'Color',colors(1,:));
% hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
plot(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),...
    'Color',colors(1,:),'DisplayName','Traj')
hpt2 = plot(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(r_encounter.EA(1),r_encounter.EA(2),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
% plot(r_encounter.astA1(1),r_encounter.astA1(2),...
%     '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
% plot(r_encounter.astD1(1),r_encounter.astD1(2),...
%     '*','Color',colors(3,:),'DisplayName',sol.asteroid_1)
% plot(r_encounter.astA2(1),r_encounter.astA2(2),...
%     '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot(r_encounter.astAa(1),r_encounter.astAa(2),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_a)
plot(r_encounter.astDa(1),r_encounter.astDa(2),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_a)
plot(r_encounter.astAb(1),r_encounter.astAb(2),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_b)

% quiver(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),2,...
%     'Color',colors(3,:),'DisplayName',sol.asteroid_1)
% quiver(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),2,...
%     'Color',colors(4,:),'DisplayName',sol.asteroid_2)
quiver(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),T_transf_orbit_a(:,1),T_transf_orbit_a(:,2),2,...
    'Color',colors(5,:),'DisplayName',sol.asteroid_a)
quiver(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),T_transf_orbit_b(:,1),T_transf_orbit_b(:,2),2,...
    'Color',colors(6,:),'DisplayName',sol.asteroid_b)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]');
% title('In-plane Thrust(shall be tangential)')
% Sun
plot(0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

Thrust_Helio_SC1 = [T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3);...
                    T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3)];
Thrust_Helio_SC2 = [T_transf_orbit_a(:,1),T_transf_orbit_a(:,2),T_transf_orbit_a(:,3);...
                    T_transf_orbit_b(:,1),T_transf_orbit_b(:,2),T_transf_orbit_b(:,3)];
                