% main 2RL GA ball simplified with rel theta
% method of JP for Dep/Lamb/LT/Lamb/DSM/Lamb/LT/Lamb/GA/LT/Ast1/LT/Ast2
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
sim.n_sol     = 100;                    % number of computational nodes
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
sim.dV_DSM_max = 0.08/sim.DU*sim.TU; % 80 m/s max dsm amount of dV allowed

%% Boundaries
% -- SC2 -- %
% dep_time fixed from SC1
% 1st segment of lambert: dep -> P1
% rP1 (1)
bound.rP1_min = 0.8; % AU, +/- 20% 1AU of Earth
bound.rP1_max = 1.2; % AU
% thetaP1 (2)
bound.thetaP1_min = 0; % rad
bound.thetaP1_max = pi; % rad, we start from a condition in the "upper" quadrant
% TOF_P1 (3), from earth dep to P1 time
bound.TOF_P1_min = 0; % days
bound.TOF_P1_max = 0.25*365*3600*24/sim.TU; % days

% 1st segment of LT: from  P1 -> P2
% Nrev fixed at 0 because i say so
% position of start is determined by end of previous leg
% position of end is once again a variable
% rP2 (4)
bound.rP2_min = 0.8; % AU, +/- 20% 1AU of Earth
bound.rP2_max = 1.2; % AU
% thetaP2 (5)
bound.thetaP2_min = 0; % rad
bound.thetaP2_max = pi/2; % rad
% TOF_P2 (6), from P1 to P2
bound.TOF_P2_min = 0*365*3600*24/sim.TU; 
bound.TOF_P2_max = 0.25*365*3600*24/sim.TU;

% 2nd segment of Lambert: from P2 -> DSM
% rDSM (7)
bound.rDSM_min = 0.8; % AU, +/- 20% 1AU of Earth
bound.rDSM_max = 1.2; % AU
% thetaDSM (8)
bound.thetaDSM_min = 0; % rad
bound.thetaDSM_max = pi/2; % rad
% TOF_DSM (9), from P2 to DSM
bound.TOF_DSM_min = 0*365*3600*24/sim.TU; 
bound.TOF_DSM_max = 0.25*365*3600*24/sim.TU;

% 3rd segment of Lambert: from DSM1 -> P3
% rP3 (10)
bound.rP3_min = 0.8; % AU, +/- 20% 1AU of Earth
bound.rP3_max = 1.2; % AU
% thetaP3 (11)
bound.thetaP3_min = 0; % rad
bound.thetaP3_max = pi/2; % rad
% TOF_P3 (12), from DSM1 to P3
bound.TOF_DSM_min = 0*365*3600*24/sim.TU; 
bound.TOF_DSM_max = 0.25*365*3600*24/sim.TU;

% 2nd segment of LT: from P3 -> P4
% Nrev fixed at 0 because i say so
% position of start is determined by end of previous leg
% position of end is once again a variable
% rP4 (13)
bound.rP4_min = 0.8; % AU, +/- 20% 1AU of Earth
bound.rP4_max = 1.2; % AU
% thetaP4 (14)
bound.thetaP4_min = 0; % rad
bound.thetaP4_max = pi/2; % rad
% TOF_P4 (15), from P3 to P4
bound.TOF_P4_min = 0*365*3600*24/sim.TU; 
bound.TOF_P4_max = 0.25*365*3600*24/sim.TU;

% 4th segment of Lambert: from P4 -> GA
% GA Stuff SC2, ballistic
% TOF_GA (16)
bound.TOFGA_a_min = 0*365*3600*24/sim.TU; 
bound.TOFGA_a_max = 1/6*365*3600*24/sim.TU; 
% azimuth out (POST-GRAVITY ASSIST)(17)
bound.az_out_GA_a_min = -pi;
bound.az_out_GA_a_max =  pi;
% elevation (POST-GRAVITY ASSIST)(18)
bound.el_out_GA_a_min = -pi/2;
bound.el_out_GA_a_max =  pi/2;

% From gravity assist SC2 to ast a, start of LT powering
% TOFa (19)
bound.TOFa_min = 0.2*365*3600*24/sim.TU; %0.2 
bound.TOFa_max = 3*365*3600*24/sim.TU; 
% N REV a (20)
bound.N_REVa_min = 0; %0
bound.N_REVa_max = 2; %3
% Coasting time a (21)
bound.CTa_min = 30*3600*24/sim.TU; % 30 days
bound.CTa_max = 80*3600*24/sim.TU; 
% From asta+CTa to ast b
% TOFb (22)
bound.TOFb_min = 0.2*365*3600*24/sim.TU; %0.2
bound.TOFb_max = 3*365*3600*24/sim.TU; 
% N REV b (23)
bound.N_REVb_min = 0; %0
bound.N_REVb_max = 2; %3

% % -- ID Permutation (24)
bound.IDP2_min = 1; 
bound.IDP2_max = data.HowMany; 

% retrieve 1st sc stuff
% Departure dates (25)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;

% out of ecliptic angle for let NLI work
% iP1 (26)
bound.iP1_min = -deg2rad(0.5); % rad
bound.iP1_max = deg2rad(0.5); % rad
% iP2 (27)
bound.iP2_min = -deg2rad(0.5); % rad
bound.iP2_max = deg2rad(0.5); % rad
% iDSM (28)
bound.iDSM_min = -deg2rad(0.5); % rad
bound.iDSM_max = deg2rad(0.5); % rad
% iP3 (29)
bound.iP3_min = -deg2rad(0.5); % rad
bound.iP3_max = deg2rad(0.5); % rad
% iP4 (30)
bound.iP4_min = -deg2rad(0.5); % rad
bound.iP4_max = deg2rad(0.5); % rad

% x = [
% rP1 (1)
% thetaP1 (2)
% TOF_P1 (3), from earth dep to P1 time
% rP2 (4)
% thetaP2 (5)
% TOF_P2 (6), from P1 to P2
% rDSM (7)
% thetaDSM (8)
% TOF_DSM (9), from P2 to DSM
% rP3 (10)
% thetaP3 (11)
% TOF_P3 (12), from DSM1 to P3
% rP4 (13)
% thetaP4 (14)
% TOF_P4 (15), from P3 to P4
% TOF_GA (16)
% azimuth out (POST-GRAVITY ASSIST)(17)
% elevation (POST-GRAVITY ASSIST)(18)
% TOFa (19)
% N REV a (20)
% Coasting time a (21)
% TOFb (22)
% N REV b (23)
% % -- ID Permutation (24)]

% for now we fix asteroids!!!
bound.lb = [bound.rP1_min,bound.thetaP1_min,bound.TOF_P1_min,bound.rP2_min,bound.thetaP2_min,...
            bound.TOF_P2_min,bound.rDSM_min,bound.thetaDSM_min,bound.TOF_DSM_min,...
            bound.rP3_min,bound.thetaP3_min,bound.TOF_DSM_min,bound.rP4_min,bound.thetaP4_min,...
            bound.TOF_P4_min,bound.TOFGA_a_min,bound.az_out_GA_a_min,bound.el_out_GA_a_min,...
            bound.TOFa_min,bound.N_REVa_min,bound.CTa_min,bound.TOFb_min,bound.N_REVb_min,...
            bound.IDP2_min,bound.mjd2000_ed,bound.iP1_min,bound.iP2_min,bound.iDSM_min,...
            bound.iP3_min,bound.iP4_min];
bound.ub = [bound.rP1_max,bound.thetaP1_max,bound.TOF_P1_max,bound.rP2_max,bound.thetaP2_max,...
            bound.TOF_P2_max,bound.rDSM_max,bound.thetaDSM_max,bound.TOF_DSM_max,...
            bound.rP3_max,bound.thetaP3_max,bound.TOF_DSM_max,bound.rP4_max,bound.thetaP4_max,...
            bound.TOF_P4_max,bound.TOFGA_a_max,bound.az_out_GA_a_max,bound.el_out_GA_a_max,...
            bound.TOFa_max,bound.N_REVa_max,bound.CTa_max,bound.TOFb_max,bound.N_REVb_max,...
            bound.IDP2_max,bound.mjd2000_ld,bound.iP1_max,bound.iP2_max,bound.iDSM_max,...
            bound.iP3_max,bound.iP4_max];
% bound.lb2 = [bound.thetaD2_min,bound.TOFGA_a_min,...
%             bound.az_out_GA_a_min,bound.el_out_GA_a_min,bound.TOFa_min,...
%             bound.N_REVa_min,bound.CTa_min,bound.TOFb_min,bound.N_REVb_min,bound.IDP2_min];
% bound.ub2 = [bound.thetaD2_max,bound.TOFGA_a_max,...
%             bound.az_out_GA_a_max,bound.el_out_GA_a_max,bound.TOFa_max,...
%             bound.N_REVa_max,bound.CTa_max,bound.TOFb_max,bound.N_REVb_max,bound.IDP2_max];

%output = interpoliamo_insieme(output)  

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (20) NREVa, (23) NREVb, (24) IDP
constr.IntCon = [20,23,24];

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
options.MaxStallGenerations = ceil(options.MaxGenerations/30);

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
% FitnessFunction2 = @(x2) ff_2RL_GA_ball_simplified_2nd_SC(x2,sim,data,power_propulsion_data); % Function handle to the fitness function
FitnessFunction = @(x) ff_2RL_GA_ball_simplified_JP_rel_theta(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Building the solution structure
sol.departure_date = mjd20002date(x(25)*sim.TU/86400);
sol.departure_mjd2000 = x(25)*sim.TU/86400;
% 1st spacecraft characteristic times
sol.TOF_P1 = x(3); % from earth dep to P1 time
sol.TOF_P2 = x(6); % from P1 to P2
sol.TOF_DSM = x(9); % from P2 to DSM
sol.TOF_P3 = x(12); % from DSM1 to P3
sol.TOF_P4 = x(15); % from P3 to P4
sol.TOF_GAa = x(16);
sol.TOFa = x(19);
sol.CTa = x(21);
sol.TOFb = x(22);

sol.rP1_magn = x(1);
sol.thetaP1 = x(2);
sol.rP2_magn = x(4);
sol.rel_thetaP2 = x(5);
sol.rDSM_magn = x(7);
sol.rel_thetaDSM = x(8);
sol.rP3_magn = x(10);
sol.rel_thetaP3 = x(11);
sol.rP4_magn = x(13);
sol.rel_thetaP4 = x(14);
sol.iP1 = x(26);
sol.iP2 = x(27);
sol.iDSM = x(28);
sol.iP3 = x(29);
sol.iP4 = x(30);

sol.az_out_GAa = x(17);
sol.el_out_GAa = x(18);

sol.N_reva = x(20);
sol.N_revb = x(23);

%% characteristic quantities plot
[output, r_encounter, v_encounter, sol] = plot_ff_2RL_GA_ball_simplified_JP_rel_theta(x,sim,data,sol);

%% only mass e thrust
subplot(2,1,1)
plot(output.t_SC2*sim.TU/86400,output.m_SC2,'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.m_SC2,'Color',colors(2,:));
xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(1,:)); 
xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
% xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(2,:)); 
% xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
% xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
xlabel('Time [days]')
ylabel('Mass [kg]')

subplot(2,1,2)
hp1 = plot(output.t_SC2*sim.TU/86400,output.T_magn_SC2,'Color',colors(1,:));
xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(1,:)); 
xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
% xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(2,:)); 
% xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
% xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
xlabel('Time [days]')
ylabel('Thrust [N]')

lgd = legend(hp1,'SC2');
% lgd.Location = 'eastoutside';
lgd.NumColumns = 1;
clearvars lgd hp1

%% orbit plots
% transfer orbits
% SC2
r_transf_orbit_a  = [output.r.lega.*cos(output.theta.lega), ...
    output.r.lega.*sin(output.theta.lega), output.z.lega];
R_transf_orbit_a = rotate_local2ecplitic(r_encounter.GAa,r_transf_orbit_a,sim.n_sol,output.Href.lega);

r_transf_orbit_b  = [output.r.legb.*cos(output.theta.legb), ...
    output.r.legb.*sin(output.theta.legb), output.z.legb];
R_transf_orbit_b = rotate_local2ecplitic(r_encounter.astDa,r_transf_orbit_b,sim.n_sol,output.Href.legb);

% SC 2-First leg: Earth -> DSM
departure_earth_sec = sol.departure_mjd2000*3600*24;
arrival_dsm2_sec = (sol.departure_mjd2000+sol.TOF_DSM2)*3600*24;
y0_EA_DSM2 = [sim.r_EA_dep*sim.DU; sim.V_EA_DSM_dep*sim.DU/sim.TU]; %km, km/s; velocity from lambert arc transfer orbit injection
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_EA_DSM2,y_EA_DSM2] = ode113(@rates, [departure_earth_sec arrival_dsm2_sec], y0_EA_DSM2,options,'sun');
% SC 1-First leg: Earth -> DSM
departure_dsm2_sec = (sol.departure_mjd2000+sol.TOF_DSM2)*3600*24;
arrival_GAa_sec = (sol.departure_mjd2000+sol.TOF_DSM2+sol.TOFGA_a)*3600*24;
y0_DSM2_GA = [r_encounter.DSM2*sim.DU; v_encounter.dep_DSM2_GAa*sim.DU/sim.TU]; %km, km/s; velocity from lambert arc transfer orbit injection
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_DSM2_GA,y_DSM2_GA] = ode113(@rates, [departure_dsm2_sec arrival_GAa_sec], y0_DSM2_GA,options,'sun');

figure()
% plot3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),...
%     'Color',colors(1,:),'DisplayName','Traj SC2')
hold on
% plot3( y_EA_DSM2(:,1)./sim.DU, y_EA_DSM2(:,2)./sim.DU, y_EA_DSM2(:,3)./sim.DU,'Color',colors(2,:),...
%     'DisplayName','EA to DSM2');
% hdsm2 = plot3(y_DSM2_GA(:,1)./sim.DU, y_DSM2_GA(:,2)./sim.DU, y_DSM2_GA(:,3)./sim.DU,...
%     'Color',colors(3,:),'DisplayName','DSM2 to GA');
hptb = plot3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),...
    'Color',colors(1,:));
hptb.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(sim.r_EA_dep(1),sim.r_EA_dep(2),sim.r_EA_dep(3),...
    'o','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.DSM2(1),r_encounter.DSM2(2),r_encounter.DSM2(3),...
    '^','Color',colors(5,:),'DisplayName','DSM')
plot3(r_encounter.GAa(1),r_encounter.GAa(2),r_encounter.GAa(3),...
    'd','Color',colors(8,:),'DisplayName','Earth GA')
plot3(r_encounter.astAa(1),r_encounter.astAa(2),r_encounter.astAa(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Arr')
plot3(r_encounter.astDa(1),r_encounter.astDa(2),r_encounter.astDa(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_a+' Dep')
plot3(r_encounter.astAb(1),r_encounter.astAb(2),r_encounter.astAb(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_b+' Arr')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(sim.departure_adim*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(sim.departure_adim*sim.TU/(3600*24),4,colors,6); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+output.ta(end)*sim.TU/(3600*24);
hello_orbit2 = sol.departure_mjd2000+(output.ta(end)+ sol.CTa + output.tb(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_a,colors,3);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_b,colors,4);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)