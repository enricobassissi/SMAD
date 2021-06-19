% main 2RL GA ball simplified part2
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

load(imp_path+'Asteroids\Workspaces\ws_2RL_0.27_GA_Ball.mat')


% % 2ND SPACECRAFT ASTEROID OBJECTIVES
asteroid_sequence = [sol.asteroid_1,sol.asteroid_2];
TF = contains(data.asteroid_names,asteroid_sequence);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, data.PermutationMatrix_2SC, data.HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);

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
sim.M2 = 100; % SC wet mass [kg] %%
sim.M_pods = 5; % mass of the payloads + landing stuff [kg] + propellant
sim.max_Available_Thrust = 0.015; % 10 [mN], BepiColombo is 250 mN but it's much bigger
sim.dV_launch_man = 0.08/sim.DU*sim.TU; % 80 m/s max manouvre at beginning 
sim.C3_max = 30/(sim.DU^2)*(sim.TU^2); % max C3 launcher, km^2/s^2, adimensionalised
sim.dV_DSM_max = 0.1*sim.TU/sim.DU; % max 200 m/s

%% Boundaries
% -- SC2 -- %
% DSM stuff SC2, ballistic
% thetaD (1)
bound.thetaD2_min = -pi/5; % rad, about +/- 10° wrt first DSM
bound.thetaD2_max = pi/5; % rad

% GA Stuff SC2, ballistic
% TOFGA (2), from DSM1 to GA 
bound.TOFGA_a_min = 0.1*365*3600*24/sim.TU; 
bound.TOFGA_a_max = 1.2*365*3600*24/sim.TU; 
% azimuth out (POST-GRAVITY ASSIST)(3)
bound.az_out_GA_a_min = -pi;
bound.az_out_GA_a_max =  pi;
% elevation (POST-GRAVITY ASSIST)(4)
bound.el_out_GA_a_min = -pi/2;
bound.el_out_GA_a_max =  pi/2;

% From gravity assist SC2 to ast a, start of LT powering
% TOFa (5)
bound.TOFa_min = 0.2*365*3600*24/sim.TU; %0.2 
bound.TOFa_max = 3*365*3600*24/sim.TU; 
% N REV a (6)
bound.N_REVa_min = 0; %0
bound.N_REVa_max = 2; %3
% Coasting time a (7)
bound.CTa_min = 30*3600*24/sim.TU; % 30 days
bound.CTa_max = 80*3600*24/sim.TU; 
% From asta+CTa to ast b
% TOFb (8)
bound.TOFb_min = 0.2*365*3600*24/sim.TU; %0.2
bound.TOFb_max = 3*365*3600*24/sim.TU; 
% N REV b (9)
bound.N_REVb_min = 0; %0
bound.N_REVb_max = 2; %3

% % -- ID Permutation (10)
bound.IDP2_min = 1; 
bound.IDP2_max = data.HowMany_2SC; 

% x = [
% thetaD2 (1)
% TOFGAa (2), from DSM2 to GA
% azimuth out (POST-GRAVITY ASSIST)(3)
% elevation (POST-GRAVITY ASSIST)(4)
% TOFa (5)
% N REV a (6)
% Coasting time a (7)
% TOFb (8)
% N REV b (9)

% -- ID Permutation2 (10)]

% for now we fix asteroids!!!
% bound.lb = [bound.mjd2000_ed,bound.TOF_DSM1_min,bound.rD1_min,bound.thetaD1_min,bound.TOFGA_1_min,...
%             bound.az_out_GA_1_min,bound.el_out_GA_1_min,bound.TOF1_min,...
%             bound.N_REV1_min,bound.CT1_min,bound.TOF2_min,bound.N_REV2_min,bound.IDP_min];
% bound.ub = [bound.mjd2000_ld,bound.TOF_DSM1_max,bound.rD1_max,bound.thetaD1_max,bound.TOFGA_1_max,...
%             bound.az_out_GA_1_max,bound.el_out_GA_1_max,bound.TOF1_max,...
%             bound.N_REV1_max,bound.CT1_max,bound.TOF2_max,bound.N_REV2_max,bound.IDP_max];
bound.lb2 = [bound.thetaD2_min,bound.TOFGA_a_min,...
            bound.az_out_GA_a_min,bound.el_out_GA_a_min,bound.TOFa_min,...
            bound.N_REVa_min,bound.CTa_min,bound.TOFb_min,bound.N_REVb_min,bound.IDP2_min];
bound.ub2 = [bound.thetaD2_max,bound.TOFGA_a_max,...
            bound.az_out_GA_a_max,bound.el_out_GA_a_max,bound.TOFa_max,...
            bound.N_REVa_max,bound.CTa_max,bound.TOFb_max,bound.N_REVb_max,bound.IDP2_max];


%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (6) NREVa, (9) NREVb, (10) IDP
constr.IntCon2 = [6,9,10];

% date of departure
sim.departure_adim = x(1);
sim.theta_DSM1 = sol.thetaD1;
sim.V_EA_DSM_dep = v_encounter.dep_EA_DSM1;
sim.r_EA_dep = r_encounter.EA;

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
FitnessFunction2 = @(x2) ff_2RL_GA_ball_simplified_2nd_SC(x2,sim,data,power_propulsion_data); % Function handle to the fitness function
numberOfVariables2 = length(bound.ub2); % Number of decision variables

tic
[x2,Fval2,exitFlag2,Output2] = ga(FitnessFunction2,numberOfVariables2,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb2,...
    bound.ub2,constr.nonlcon,constr.IntCon2,options);
el_time_min_pp2 = toc/60;

%% Building the solution structure
sol2.departure_date = mjd20002date(sim.departure_adim*sim.TU/86400);
sol2.departure_mjd2000 = sim.departure_adim*sim.TU/86400;
% 1st spacecraft characteristic times
sol2.TOFGA_a = x2(2)*sim.TU/86400; % DSM to GA
sol2.TOFa = x2(5)*sim.TU/86400; % tof sc1 GA to Ast1 
sol2.CTa = x2(7)*sim.TU/86400;
sol2.TOFb = x2(8)*sim.TU/86400; % tof sc1 from ast 1 to 2nd asteroid
% mission duration
% sol2.mission_duration_d_SC2 = sol2.TOF_DSM2+sol2.TOFGA_a+sol2.TOFa+sol2.CT1+sol2.TOFb;
% sol2.mission_duration_y_SC2 = sol2.mission_duration_d_SC2/365;
% sol2.arrival_dat = mjd20002date(sim.departure_adim*sim.TU/(3600*24)+sol2.mission_duration_d_SC2);
% -- DSM SC2
sol2.thetaD2 = x2(1); % rad
% asteorids
sol2.IDP2 = x2(10);
sol2.asteroid_a = data.PermutationMatrix_2SC(sol2.IDP2,1);
sol2.asteroid_b = data.PermutationMatrix_2SC(sol2.IDP2,2);
% sol2.asteroid_a = "2010UJ";
% sol2.asteroid_b = "2016CK137";

sol2.Nrev = [x2(6), x2(9)];
% GA stuff SC2 - OUT
sol2.az_out_GAa = x(3);
sol2.el_out_GAa = x(4);

%% characteristic quantities plot
[output2, r_encounter2, v_encounter2, sol2] = plot_ff_2RL_GA_ball_simplified_2nd_SC(x2,sim,data,sol2,power_propulsion_data);

%% only mass e thrust
subplot(2,1,1)
plot(output2.t_SC2*sim.TU/86400,output2.m_SC2,'Color',colors(1,:));
% hold on
% plot(output.t_SC2*sim.TU/86400,output.m_SC2,'Color',colors(2,:));
xline(sol2.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(1,:)); 
xline(sol2.TOFa+sol2.CTa,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
xline(sol2.TOFa+sol2.CTa+sol2.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
% xline(sol.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(2,:)); 
% xline(sol.TOFa+sol.CTa,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
% xline(sol.TOFa+sol.CTa+sol.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(2,:));
xlabel('Time [days]')
ylabel('Mass [kg]')

subplot(2,1,2)
hp1 = plot(output2.t_SC2*sim.TU/86400,output2.T_magn_SC2,'Color',colors(1,:));
xline(sol2.TOFa,'LineWidth',2,'LineStyle',':','Color',colors(1,:)); 
xline(sol2.TOFa+sol2.CTa,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
xline(sol2.TOFa+sol2.CTa+sol2.TOFb,'LineWidth',2,'LineStyle',':','Color',colors(1,:));
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
r_transf_orbit_a  = [output2.r.lega.*cos(output2.theta.lega), ...
    output2.r.lega.*sin(output2.theta.lega), output2.z.lega];
R_transf_orbit_a = rotate_local2ecplitic(r_encounter2.GAa,r_transf_orbit_a,sim.n_sol,output2.Href.lega);

r_transf_orbit_b  = [output2.r.legb.*cos(output2.theta.legb), ...
    output2.r.legb.*sin(output2.theta.legb), output2.z.legb];
R_transf_orbit_b = rotate_local2ecplitic(r_encounter2.astDa,r_transf_orbit_b,sim.n_sol,output2.Href.legb);

% SC 2-First leg: Earth -> DSM
departure_earth_sec = sol2.departure_mjd2000*3600*24;
arrival_dsm2_sec = (sol2.departure_mjd2000+sol2.TOF_DSM2)*3600*24;
y0_EA_DSM2 = [sim.r_EA_dep*sim.DU; sim.V_EA_DSM_dep*sim.DU/sim.TU]; %km, km/s; velocity from lambert arc transfer orbit injection
options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_EA_DSM2,y_EA_DSM2] = ode113(@rates, [departure_earth_sec arrival_dsm2_sec], y0_EA_DSM2,options,'sun');
% SC 1-First leg: Earth -> DSM
departure_dsm2_sec = (sol2.departure_mjd2000+sol2.TOF_DSM2)*3600*24;
arrival_GAa_sec = (sol2.departure_mjd2000+sol2.TOF_DSM2+sol2.TOFGA_a)*3600*24;
y0_DSM2_GA = [r_encounter2.DSM2*sim.DU; v_encounter2.dep_DSM2_GAa*sim.DU/sim.TU]; %km, km/s; velocity from lambert arc transfer orbit injection
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
plot3(r_encounter2.DSM2(1),r_encounter2.DSM2(2),r_encounter2.DSM2(3),...
    '^','Color',colors(5,:),'DisplayName','DSM')
plot3(r_encounter2.GAa(1),r_encounter2.GAa(2),r_encounter2.GAa(3),...
    'd','Color',colors(8,:),'DisplayName','Earth GA')
plot3(r_encounter2.astAa(1),r_encounter2.astAa(2),r_encounter2.astAa(3),...
    '^','Color',colors(3,:),'DisplayName',sol2.asteroid_a+' Arr')
plot3(r_encounter2.astDa(1),r_encounter2.astDa(2),r_encounter2.astDa(3),...
    '*','Color',colors(3,:),'DisplayName',sol2.asteroid_a+' Dep')
plot3(r_encounter2.astAb(1),r_encounter2.astAb(2),r_encounter2.astAb(3),...
    '^','Color',colors(4,:),'DisplayName',sol2.asteroid_b+' Arr')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(sim.departure_adim*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(sim.departure_adim*sim.TU/(3600*24),4,colors,6); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol2.departure_mjd2000+output2.ta(end)*sim.TU/(3600*24);
hello_orbit2 = sol2.departure_mjd2000+(output2.ta(end)+ sol2.CTa + output2.tb(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol2.asteroid_a,colors,3);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol2.asteroid_b,colors,4);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)