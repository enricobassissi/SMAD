%% --------------------------------------------------------------------- %%
%% --------------- Earth Ast1 Ast2 Ast3 Ast4 Transfer ------------------ %%
%% ------------------------- ARCH 1FL ---------------------------------- %%
%% --------------------- FREEDOM AT EACH AST --------------------------- %%
%% ------------------------------ SOO ---------------------------------- %%
%% --------------------------------------------------------------------- %%
%% Setup for default options
set(0, 'DefaultTextFontSize',11)
set(0, 'DefaultAxesFontSize', 11)
set(0, 'DefaultLegendFontSize', 11)
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.2)
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

sim.case_name = 'ARCH 2FL: 2SC FB LOW THRUST EACH ON TWO ASTEROID';

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

load('data_processed_42_61_2SC.mat')


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
sim.M1 = 70; % SC wet mass [kg]
sim.M2 = 70; % SC wet mass [kg]
sim.M_pods = 5; % mass of the payloads+landing stuff [kg]

sim.c = 1e+3; %case 12: 1e+3
sim.max_Available_Thrust = 0.01; %case 12: 0.02

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; %2024
bound.date_ld = [2028, 1, 1, 0, 0, 0]; %2028
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 0.3*3600*24/sim.TU; %0.3*365
bound.TOF1_max = 3*365*3600*24/sim.TU; %2*365
% TOF2 (3)
bound.TOF2_min = 0.3*3600*24/sim.TU; 
bound.TOF2_max = 3*365*3600*24/sim.TU; 
% TOFa (4)
bound.TOFa_min = 0.3*365*3600*24/sim.TU; 
bound.TOFa_max = 3*365*3600*24/sim.TU; 
% TOFb (5)
bound.TOFb_min = 0.3*365*3600*24/sim.TU; 
bound.TOFb_max = 3*365*3600*24/sim.TU; 
% N REV 1 (6)
bound.N_REV1_min = 0; %0
bound.N_REV1_max = 2; %3
% N REV 2 (7)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 2; %3
% N REVa (8)
bound.N_REVa_min = 0; %0
bound.N_REVa_max = 2; %3
% N REVb (9)
bound.N_REVb_min = 0; %0
bound.N_REVb_max = 2; %3
% ID Permutation (10)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 
% ID Permutation 2 (11)
% it's 1-100 because then inside it gets rescaled depending on the HowMany
% valid sequence of asteroids are left depending on the first 2 chosen
bound.IDP2_min = 1; 
bound.IDP2_max = 100; 

% C3 stuff
% Constraint on C3 Launcher (12)
sim.C3_max = 30; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (13)
bound.az_min = -pi;
bound.az_max = pi;
% elevation (14)
bound.el_min = -pi/2;
bound.el_max = pi/2;

% FB on asteroids
% ASTEROID 1 FLYBY
% vinf_ast1 (15)
bound.vinf_ast1_min = 0;
bound.vinf_ast1_max = 2*sim.TU/sim.DU;
% azimuth_ast1_in (16)
bound.az_ast1_min = -pi;
bound.az_ast1_max = pi;
% elevation_ast1_in (17)
bound.el_ast1_min = -pi/2;
bound.el_ast1_max = pi/2;

% ASTEROID 2 FLYBY
% vinf_ast2 (18)
bound.vinf_ast2_min = 0;
bound.vinf_ast2_max = 2*sim.TU/sim.DU;
% azimuth_ast2_in (19)
bound.az_ast2_min = -pi;
bound.az_ast2_max = pi;
% elevation_ast2_in (20)
bound.el_ast2_min = -pi/2;
bound.el_ast2_max = pi/2;

% ASTEROID a FLYBY
% vinf_asta (21)
bound.vinf_asta_min = 0;
bound.vinf_asta_max = 2*sim.TU/sim.DU;
% azimuth_asta_in (22)
bound.az_asta_min = -pi;
bound.az_asta_max = pi;
% elevation_asta_in (23)
bound.el_asta_min = -pi/2;
bound.el_asta_max = pi/2;

% ASTEROID b FLYBY
% vinf_astb (24)
bound.vinf_astb_min = 0;
bound.vinf_astb_max = 2*sim.TU/sim.DU;
% azimuth_astb_in (25)
bound.az_astb_min = -pi;
bound.az_astb_max = pi;
% elevation_astb_in (26)
bound.el_astb_min = -pi/2;
bound.el_astb_max = pi/2;

% x = [MJD0,TOF1,TOF2,TOFa,TOFb,NREV1,NREV2,NREVa,NREVb,IDP, IDP2,v_inf_magn,az,el,] and all
%      asteroids flyby variables

bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.TOFa_min, bound.TOFb_min, bound.N_REV1_min, bound.N_REV2_min, ...
    bound.N_REVa_min, bound.N_REVb_min, bound.IDP_min,bound.IDP2_min, ...
    bound.v_inf_magn_min, bound.az_min, bound.el_min,  ...
    bound.vinf_ast1_min, bound.az_ast1_min, bound.el_ast1_min, ...
    bound.vinf_ast2_min, bound.az_ast2_min, bound.el_ast2_min, ...
    bound.vinf_asta_min, bound.az_asta_min, bound.el_asta_min, ...
    bound.vinf_astb_min, bound.az_astb_min, bound.el_astb_min]; % Lower bound


bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.TOFa_max, bound.TOFb_max, bound.N_REV1_max, bound.N_REV2_max, ...
    bound.N_REVa_max, bound.N_REVb_max, bound.IDP_max, bound.IDP2_min, ...
    bound.v_inf_magn_max, bound.az_max, bound.el_max,  ...
    bound.vinf_ast1_max, bound.az_ast1_max, bound.el_ast1_max, ...
    bound.vinf_ast2_max, bound.az_ast2_max, bound.el_ast2_max, ...
    bound.vinf_asta_max, bound.az_asta_max, bound.el_asta_max, ...
    bound.vinf_astb_max, bound.az_astb_max, bound.el_astb_max]; % Upper bound


%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (6) NREV1, (7) NREV2, (8) NREV3, (9) NREV4, (10) IDP, (11) IDP2
constr.IntCon = [6,7,8,9,10,11];

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
FitnessFunction = @(x) ff_2FL(x,sim,data); % Function handle to the fitness function
numberOfVariables = length(bound.ub); % Number of decision variables

tic
[x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,constr.A, ...
    constr.b,constr.Aeq,constr.beq,bound.lb,...
    bound.ub,constr.nonlcon,constr.IntCon,options);
el_time_min_pp = toc/60;

%% Building the solution structure
sol.asteroid_1 = data.PermutationMatrix(x(10),1);
sol.asteroid_2 = data.PermutationMatrix(x(10),2);

IDP_temp_2 = x(11); % index for 2nd permutation matrix to be built inside depending on the first 2 selected asteroids
TF = contains(data.asteroid_names,[sol.asteroid_1,sol.asteroid_1]);
data_elements_matrix_2SC = data.data_elements_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning2(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
sol.asteroid_a = PermutationMatrix_2SC(IDP2,1);
sol.asteroid_b = PermutationMatrix_2SC(IDP2,2);

sol.departure_mjd2000 = x(1)*sim.TU/(3600*24);
sol.departure_date = mjd20002date(sol.departure_mjd2000);
sol.TOF1  = x(2)*sim.TU/(3600*24);
sol.TOF2  = x(3)*sim.TU/(3600*24);
sol.TOFa  = x(4)*sim.TU/(3600*24);
sol.TOFb  = x(5)*sim.TU/(3600*24);
sol.mission_duration_d_SC1 = sol.TOF1+sol.TOF2;
sol.mission_duration_d_SC2 = sol.TOFa+sol.TOFb;
sol.mission_duration_y_SC1 = sol.mission_duration_d_SC1/365;
sol.mission_duration_y_SC2 = sol.mission_duration_d_SC2/365;

sol.arrival_date_latest = mjd20002date(x(1)*sim.TU/(3600*24)... 
                          + max(sol.mission_duration_d_SC1,sol.mission_duration_d_SC2));

sol.objfun = Fval;
sol.Nrev = [x(6), x(7), x(8), x(9)];
sol.v_inf_magn = x(12)*sim.DU/sim.TU;
sol.az = x(13);
sol.el = x(14);

sol.v_inf_ast1_magn = x(15)*sim.DU/sim.TU;
sol.v_inf_ast2_magn = x(18)*sim.DU/sim.TU;
sol.v_inf_asta_magn = x(21)*sim.DU/sim.TU;
sol.v_inf_astb_magn = x(24)*sim.DU/sim.TU;

%% characteristic quantities plot and Checks
[output, r_encounter, sol] = plot_ff_2FL(x,sim,data,sol);

figure()
subplot(5,1,1)
hp1 = plot(output.t_SC1*sim.TU/86400,output.Thrust_SC1(:,1),'Color',colors(1,:));
hold on
hp2 = plot(output.t_SC2*sim.TU/86400,output.Thrust_SC2(:,1),'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xline(output.t(sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(2*sim.n_sol)*sim.TU/86400,'LineWidth',2);
% xline(output.t(3*sim.n_sol)*sim.TU/86400,'LineWidth',2); xline(output.t(4*sim.n_sol)*sim.TU/86400,'LineWidth',2);
% xline(output.t(5*sim.n_sol)*sim.TU/86400,'LineWidth',2)
% xlabel('Time [days]')
ylabel('In-plane Thrust [N]')
lgd = legend([hp1,hp2],'SC1','SC2');
lgd.Location = 'eastoutside';
lgd.NumColumns = 1;
clearvars lgd hp1 hp2

subplot(5,1,2)
plot(output.t_SC1*sim.TU/86400,output.Thrust_SC1(:,2),'Color',colors(1,:));
hold on
plot(output.t_SC2*sim.TU/86400,output.Thrust_SC2(:,2),'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('In-plane Thrust angle [deg]')

subplot(5,1,3)
plot(output.t_SC1*sim.TU/86400,output.Thrust_SC1(:,3),'Color',colors(1,:));
hold on
plot(output.t_SC2*sim.TU/86400,output.Thrust_SC2(:,3),'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('out-of-plane Thrust [N]')

subplot(5,1,4)
plot(output.t_SC1*sim.TU/86400,output.T_magn_SC1,'Color',colors(1,:));
hold on
plot(output.t_SC2*sim.TU/86400,output.T_magn_SC2,'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
% xlabel('Time [days]')
ylabel('Thrust [N]')

subplot(5,1,5)
plot(output.t_SC1*sim.TU/86400,output.m_SC1,'Color',colors(1,:));
hold on
plot(output.t_SC2*sim.TU/86400,output.m_SC2,'Color',colors(2,:));
% xline(sol.TOF1,'LineWidth',2); xline(sol.TOF1+sol.CT1,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2,'LineWidth',2);xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2,'LineWidth',2);
% xline(sol.TOF1+sol.CT1+sol.TOF2+sol.CT2+sol.TOF3,'LineWidth',2);
xlabel('Time [days]')
ylabel('Mass [kg]')

%% orbit plots
% transfer orbits
r_transf_orbit_1  = [output.r.leg1.*cos(output.theta.leg1), ...
    output.r.leg1.*sin(output.theta.leg1), output.z.leg1];
R_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_1,sim.n_sol,output.Href.leg1);

r_transf_orbit_2  = [output.r.leg2.*cos(output.theta.leg2), ...
    output.r.leg2.*sin(output.theta.leg2), output.z.leg2];
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.ast1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

r_transf_orbit_a  = [output.r.lega.*cos(output.theta.lega), ...
    output.r.lega.*sin(output.theta.lega), output.z.lega];
R_transf_orbit_a = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_a,sim.n_sol,output.Href.lega);

r_transf_orbit_b  = [output.r.legb.*cos(output.theta.legb), ...
    output.r.legb.*sin(output.theta.legb), output.z.legb];
R_transf_orbit_b = rotate_local2ecplitic(r_encounter.asta,r_transf_orbit_b,sim.n_sol,output.Href.legb);

figure()
hpt1 = plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Trajectory');
% hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt3 = plot3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),...
    'Color',colors(1,:));
hpt3.Annotation.LegendInformation.IconDisplayStyle = 'off';
hpt4 = plot3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),...
    'Color',colors(1,:));
hpt4.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(r_encounter.EA(1),r_encounter.EA(2),r_encounter.EA(3),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.ast2(1),r_encounter.ast2(2),r_encounter.ast2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.asta(1),r_encounter.asta(2),r_encounter.asta(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_a)
plot3(r_encounter.astb(1),r_encounter.astb(2),r_encounter.astb(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_b)
% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,2); % mars
% Asteroids
fraction_of_the_orbit = 1;
hello_orbit1 = sol.departure_mjd2000+output.t1(end)*sim.TU/(3600*24);
hello_orbit2 = sol.departure_mjd2000+(output.t1(end)+output.t2(end))*sim.TU/(3600*24);
hello_orbita = sol.departure_mjd2000+output.ta(end)*sim.TU/(3600*24);
hello_orbitb = sol.departure_mjd2000+(output.ta(end)+output.tb(end))*sim.TU/(3600*24);
plot_asteorid_orbit(hello_orbit1,fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(hello_orbit2,fraction_of_the_orbit,sol.asteroid_2,colors,4);
plot_asteorid_orbit(hello_orbita,fraction_of_the_orbit,sol.asteroid_a,colors,5);
plot_asteorid_orbit(hello_orbitb,fraction_of_the_orbit,sol.asteroid_b,colors,6);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
legend('show')
view(2)

%% post analysis and checks
max_ThrustSC1 = max(output.T_magn_SC1)*1000; %[mN]
max_ThrustSC2 = max(output.T_magn_SC2)*1000; %[mN]

transf1_check = abs(sol.TOF1 - output.t1(end))*sim.TU/(3600*24);
transf2_check = abs(sol.TOF2 - output.t2(end))*sim.TU/(3600*24);
transfa_check = abs(sol.TOFa - output.ta(end))*sim.TU/(3600*24);
transfb_check = abs(sol.TOFb - output.tb(end))*sim.TU/(3600*24);

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

% Propulsion for NIKITA
Propulsion.Magnitude_Thrust = [output.T_magn_SC1 output.T_magn_SC2];
Propulsion.InPlane_Thrust   = [output.Thrust_SC1(:,1) output.Thrust_SC2(:,1)];
Propulsion.gamma_angle      = [output.Thrust_SC1(:,2) output.Thrust_SC2(:,2)];
Propulsion.OutofPlane_Thrust = [output.Thrust_SC1(:,3) output.Thrust_SC2(:,3)];

%% Plot with thrust vectors
Tlocal_transf_orbit_1  = [-sol.T_1(:,1).*sin(output.theta.leg1), ...
    sol.T_1(:,1).*cos(output.theta.leg1), sol.T_1(:,3)];
T_transf_orbit_1 = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit_1,sim.n_sol,output.Href.leg1);

Tlocal_transf_orbit_2  = [-sol.T_2(:,1).*sin(output.theta.leg2), ...
    sol.T_2(:,1).*cos(output.theta.leg2), sol.T_2(:,3)];
T_transf_orbit_2 = rotate_local2ecplitic(r_encounter.ast1,Tlocal_transf_orbit_2,sim.n_sol,output.Href.leg2);

Tlocal_transf_orbit_a  = [-sol.T_a(:,1).*sin(output.theta.lega), ...
    sol.T_a(:,1).*cos(output.theta.lega), sol.T_a(:,3)];
T_transf_orbit_a = rotate_local2ecplitic(r_encounter.EA,Tlocal_transf_orbit_a,sim.n_sol,output.Href.lega);

Tlocal_transf_orbit_b  = [-sol.T_b(:,1).*sin(output.theta.legb), ...
    sol.T_b(:,1).*cos(output.theta.legb), sol.T_b(:,3)];
T_transf_orbit_b = rotate_local2ecplitic(r_encounter.asta,Tlocal_transf_orbit_b,sim.n_sol,output.Href.legb);

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
plot3(r_encounter.ast1(1),r_encounter.ast1(2),r_encounter.ast1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot3(r_encounter.ast2(1),r_encounter.ast2(2),r_encounter.ast2(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_2)
plot3(r_encounter.asta(1),r_encounter.asta(2),r_encounter.asta(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_a)
plot3(r_encounter.astb(1),r_encounter.astb(2),r_encounter.astb(3),...
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

figure('Name','Plane thrust')
hold on
% plot(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),...
%     'Color',colors(1,:),'DisplayName','Traj')
% hpt1 = plot(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),...
%     'Color',colors(1,:));
% hpt1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
plot(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),...
    'Color',colors(1,:),'DisplayName','Traj')
% hpt2 = plot(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),...
%     'Color',colors(1,:));
% hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(r_encounter.EA(1),r_encounter.EA(2),...
    '*','Color',colors(8,:),'DisplayName','Dep Earth')
plot(r_encounter.ast1(1),r_encounter.ast1(2),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1)
plot(r_encounter.ast2(1),r_encounter.ast2(2),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2)
plot(r_encounter.asta(1),r_encounter.asta(2),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_a)
plot(r_encounter.astb(1),r_encounter.astb(2),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_b)

% quiver(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),2,...
%     'Color',colors(3,:),'DisplayName',sol.asteroid_1)
% quiver(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),2,...
%     'Color',colors(4,:),'DisplayName',sol.asteroid_2)
quiver(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),T_transf_orbit_a(:,1),T_transf_orbit_a(:,2),2,...
    'Color',colors(5,:),'DisplayName',sol.asteroid_a)
% quiver(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),T_transf_orbit_b(:,1),T_transf_orbit_b(:,2),2,...
%     'Color',colors(6,:),'DisplayName',sol.asteroid_b)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]');
title('In-plane Thrust(shall be tangential)')
% Sun
plot(0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

Thrust_Helio_SC1 = [T_transf_orbit_1(:,1),T_transf_orbit_1(:,2),T_transf_orbit_1(:,3);...
                    T_transf_orbit_2(:,1),T_transf_orbit_2(:,2),T_transf_orbit_2(:,3)];
Thrust_Helio_SC2 = [T_transf_orbit_a(:,1),T_transf_orbit_a(:,2),T_transf_orbit_a(:,3);...
                    T_transf_orbit_b(:,1),T_transf_orbit_b(:,2),T_transf_orbit_b(:,3)];