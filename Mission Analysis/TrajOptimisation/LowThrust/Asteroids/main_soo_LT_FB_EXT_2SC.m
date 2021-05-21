%% --------------------------------------------------------------------- %%
%% ------------------------ Earth Ast1 Ast2  --------------------------- %%
%% ------------------------- Earth Asta Astb --------------------------- %%
%% ----------------------- ARCH 2SC, LT FB EXT ------------------------- %%
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
% load('data_elements_matrix_9.mat')
% data.p_number = 2;
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix, data.p_number);
% 
% data.data_element_matrix = data_elements_matrix;
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

load('data_processed_9_2.mat')

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
sim.M1 = 50; % SC wet mass [kg]
sim.M2 = 50; % SC wet mass [kg]
sim.M_pods = 5; % mass of the payloads+landing stuff [kg]
sim.max_Available_Thrust = 0.05; % 50 [mN], BepiColombo is 250 mN but it's much bigger

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 500*3600*24/sim.TU; %600
bound.TOF1_max = 900*3600*24/sim.TU; 
% TOF2 (3)
bound.TOF2_min = 0.5*365*3600*24/sim.TU; %600
bound.TOF2_max = 4*365*3600*24/sim.TU; 
% TOFa (4)
bound.TOFa_min = 500*3600*24/sim.TU; %600
bound.TOFa_max = 900*3600*24/sim.TU; 
% TOFb (5)
bound.TOFb_min = 0.5*365*3600*24/sim.TU; %600
bound.TOFb_max = 4*365*3600*24/sim.TU; 
% N REV 1 (6)
bound.N_REV1_min = 1; %0
bound.N_REV1_max = 1; %3
% N REV 2 (7)
bound.N_REV2_min = 0; %0
bound.N_REV2_max = 1; %3
% N REV a (8)
bound.N_REVa_min = 1; %0
bound.N_REVa_max = 1; %3
% N REV b (9)
bound.N_REVb_min = 0; %0
bound.N_REVb_max = 1; %3
% ID Permutation (10)
bound.IDP_min = 1; 
bound.IDP_max = data.HowMany; 
% ID Permutation 2 (11)
% it's 1-100 because then inside it gets rescaled depending on the HowMany
% valid sequence of asteroids are left depending on the first 2 chosen
bound.IDP2_min = 1; 
bound.IDP2_max = 100; 
% bound.IDP_max = length(data.asteroid_names);
% C3 stuff
% Constraint on C3 Launcher (12)
sim.C3_max = 20; % km^2/s^2
bound.v_inf_magn_min = 0;
bound.v_inf_magn_max = sqrt(sim.C3_max)/sim.DU*sim.TU;
% azimuth (13)
bound.az_min = -deg2rad(180);
bound.az_max = deg2rad(180);
% elevation (14)
bound.el_min = -deg2rad(90);
bound.el_max = deg2rad(90);

% FB on asteroid angles 
% ASTEROID 1 IN
% azimuth_ast1_in (15)
bound.az_ast1_in_min = -pi;
bound.az_ast1_in_max = pi;
% elevation_ast1_in (16)
bound.el_ast1_in_min = -pi/2;
bound.el_ast1_in_max = pi/2;
% ASTEROID 1 OUT
% azimuth_ast1_out (17)
bound.az_ast1_out_min = -pi;
bound.az_ast1_out_max = pi;
% elevation_ast1_out (18)
bound.el_ast1_out_min = -pi/2;
bound.el_ast1_out_max = pi/2;
% ASTEROID 2 IN
% azimuth_ast2_in (19)
bound.az_ast2_in_min = -pi;
bound.az_ast2_in_max = pi;
% elevation_ast2_in (20)
bound.el_ast2_in_min = -pi/2;
bound.el_ast2_in_max = pi/2;

% ASTEROID a IN
% azimuth_asta_in (21)
bound.az_asta_in_min = -pi;
bound.az_asta_in_max = pi;
% elevation_asta_in (22)
bound.el_asta_in_min = -pi/2;
bound.el_asta_in_max = pi/2;
% ASTEROID a OUT
% azimuth_asta_out (23)
bound.az_asta_out_min = -pi;
bound.az_asta_out_max = pi;
% elevation_asta_out (24)
bound.el_asta_out_min = -pi/2;
bound.el_asta_out_max = pi/2;
% ASTEROID b IN
% azimuth_astb_in (25)
bound.az_astb_in_min = -pi;
bound.az_astb_in_max = pi;
% elevation_astb_in (26)
bound.el_astb_in_min = -pi/2;
bound.el_astb_in_max = pi/2;

% x =
% [MJD0,TOF1,TOF2,TOF3,TOF4,NREV,NREV2,NREV3,NREV4,IDP,v_inf_magn,az,el]+FB EXT stuff
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.TOFa_min, bound.TOFb_min, bound.N_REV1_min, bound.N_REV2_min, ...
    bound.N_REVa_min, bound.N_REVb_min, bound.IDP_min, bound.IDP2_min, ...
    bound.v_inf_magn_min, bound.az_min, bound.el_min, ...
    bound.az_ast1_in_min, bound.el_ast1_in_min, bound.az_ast1_out_min, bound.el_ast1_out_min, ...
    bound.az_ast2_in_min, bound.el_ast2_in_min, ...
    bound.az_asta_in_min, bound.el_asta_in_min, bound.az_asta_out_min, bound.el_asta_out_min, ...
    bound.az_astb_in_min, bound.el_astb_in_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.TOFa_max, bound.TOFb_max, bound.N_REV1_max, bound.N_REV2_max, ...
    bound.N_REVa_max, bound.N_REVb_max, bound.IDP_max, bound.IDP2_max, ...
    bound.v_inf_magn_max, bound.az_max, bound.el_max, ...
    bound.az_ast1_in_max, bound.el_ast1_in_max, bound.az_ast1_out_max, bound.el_ast1_out_max, ...
    bound.az_ast2_in_max, bound.el_ast2_in_max, ...
    bound.az_asta_in_max, bound.el_asta_in_max, bound.az_asta_out_max, bound.el_asta_out_max, ...
    bound.az_astb_in_max, bound.el_astb_in_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (6) NREV1, (7) NREV2, (8) NREVa, (9) NREVb, (10) IDP, (11) IDP2
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
options.MaxGenerations = 50; % ideal 100

options.FunctionTolerance = 1e-6; %1e-9
options.MaxStallGenerations = ceil(options.MaxGenerations/5);

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
FitnessFunction = @(x) ff_ea_2SC_LT_FB_EXT_soo_NLI(x,sim,data); % Function handle to the fitness function
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
data_elements_matrix_2SC = data.data_element_matrix(~TF,:);
[~, PermutationMatrix_2SC, HowMany_2SC] = ...
            sequences_local_pruning(data_elements_matrix_2SC, data.p_number);
IDP2 = ceil(IDP_temp_2*HowMany_2SC/100);
sol.asteroid_a = PermutationMatrix_2SC(IDP2,1);
sol.asteroid_b = PermutationMatrix_2SC(IDP2,2);

sol.departure_date = mjd20002date(x(1)*sim.TU/(3600*24));
sol.TOF1 = x(2)*sim.TU/86400;
sol.TOF2 = x(3)*sim.TU/86400;
sol.TOFa = x(4)*sim.TU/86400;
sol.TOFb = x(5)*sim.TU/86400;
sol.TOF1_ADIM = x(2);
sol.TOF2_ADIM = x(3);
sol.TOFa_ADIM = x(4);
sol.TOFb_ADIM = x(5);
sol.mission_duration_d_SC1 = sol.TOF1+sol.TOF2;
sol.mission_duration_d_SC2 = sol.TOFa+sol.TOFb;
sol.mission_duration_y_SC1 = sol.mission_duration_d_SC1/365;
sol.mission_duration_y_SC2 = sol.mission_duration_d_SC2/365;
sol.arrival_date_latest = mjd20002date(x(1)*sim.TU/(3600*24)+max(sol.mission_duration_d_SC1,sol.mission_duration_d_SC2));

% sol.mass_fraction = Fval;
sol.Nrev = [x(6), x(7), x(8), x(9)];
sol.v_inf_magn = x(12)*sim.DU/sim.TU;
sol.az = x(13);
sol.el = x(14);

%% characteristic quantities plot
[output, r_encounter, v_encounter, sol] = plot_ff_ea_2SC_LT_FB_soo_NLI(x,sim,data,sol);

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
% xlabel('Time [days]')
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
plot3(R_transf_orbit_1(:,1),R_transf_orbit_1(:,2),R_transf_orbit_1(:,3),...
    'Color',colors(1,:),'DisplayName','Traj SC1')
hold on
hpt2 = plot3(R_transf_orbit_2(:,1),R_transf_orbit_2(:,2),R_transf_orbit_2(:,3),...
    'Color',colors(1,:));
hpt2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot3(R_transf_orbit_a(:,1),R_transf_orbit_a(:,2),R_transf_orbit_a(:,3),...
    'Color',colors(2,:),'DisplayName','Traj SC2');
hpt4 = plot3(R_transf_orbit_b(:,1),R_transf_orbit_b(:,2),R_transf_orbit_b(:,3),...
    'Color',colors(2,:));
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
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
% PLANETS
plot_planet_orbit(x(1)*sim.TU/(3600*24),3,colors,8); % earth
plot_planet_orbit(x(1)*sim.TU/(3600*24),4,colors,6); % mars
% Asteroids
fraction_of_the_orbit = 1/6;
plot_asteorid_orbit(output.t1(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_1,colors,3);
plot_asteorid_orbit(output.t2(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_2,colors,4);
plot_asteorid_orbit(output.ta(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_a,colors,5);
plot_asteorid_orbit(output.tb(end)*sim.TU/(3600*24),fraction_of_the_orbit,sol.asteroid_b,colors,6);
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
legend('show')
view(2)

