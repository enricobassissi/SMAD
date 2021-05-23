%% --------------------------------------------------------------------- %%
%% ------------------------ Earth Ast1 Ast2  --------------------------- %%
%% ------------------------- Earth Asta Astb --------------------------- %%
%% ------------------------- ARCH 2SC, LT RV --------------------------- %%
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
% load('data_elements_matrix_42_61_2SC.mat')
% data.p_number = 2;
% [data.asteroid_names, data.PermutationMatrix, data.HowMany] = ...
%             sequences_local_pruning(data_elements_matrix, data.p_number);
% 
% data.data_element_matrix = data_elements_matrix;
% 
% [data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

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
sim.max_Available_Thrust = 0.01; % 50 [mN], BepiColombo is 250 mN but it's much bigger
sim.ID_FLYBY = 3; % flyby planet

%% Boundaries
% Departure dates (1)
bound.date_ed = [2024, 1, 1, 0, 0, 0]; 
bound.date_ld =  [2028, 1, 1, 0, 0, 0]; 
bound.mjd2000_ed = date2mjd2000(bound.date_ed)*3600*24/sim.TU;
bound.mjd2000_ld = date2mjd2000(bound.date_ld)*3600*24/sim.TU;
% TOF1 (2)
bound.TOF1_min = 1*365*3600*24/sim.TU; %600
bound.TOF1_max = 3*365*3600*24/sim.TU; 
% TOF2 (3)
bound.TOF2_min = 0.5*365*3600*24/sim.TU; %600
bound.TOF2_max = 3*365*3600*24/sim.TU; 
% TOFa (4)
bound.TOFa_min = 1*365*3600*24/sim.TU; %600
bound.TOFa_max = 3*365*3600*24/sim.TU; 
% TOFb (5)
bound.TOFb_min = 0.5*365*3600*24/sim.TU; %600
bound.TOFb_max = 3*365*3600*24/sim.TU; 
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

% Coasting time 1 (15)
bound.CT1_min = 30*3600*24/sim.TU; % 30 days
bound.CT1_max = 80*3600*24/sim.TU; 
% Coasting time 2 (16)
bound.CTa_min = 30*3600*24/sim.TU; % 30 days
bound.CTa_max = 80*3600*24/sim.TU; 

% GA Stuff SC1
% TOFGA (17)
bound.TOFGA_min = 50*3600*24/sim.TU; %400
bound.TOFGA_max = 900*3600*24/sim.TU; 
% N REV GA (18)
bound.N_REVGA_min = 0; %0
bound.N_REVGA_max = 1; %3
% Incoming Leg
% v_inf_magn GA (19)
bound.v_inf_magnGA_min = 0;
bound.v_inf_magnGA_max = 15/sim.DU*sim.TU; % Longusky
% azimuth GAin (20)
bound.azGAin_min = -pi;
bound.azGAin_max = pi;
% elevation GAin (21)
bound.elGAin_min = -pi/2;
bound.elGAin_max = pi/2;
% Outgoing leg
% azimuth GAout (22)
bound.azGAout_min = -pi;
bound.azGAout_max = pi;
% elevation GAout (23)
bound.elGAout_min = -pi/2;
bound.elGAout_max = pi/2;

% GA Stuff SC2
% TOFGA (24)
bound.TOFGA2_min = 50*3600*24/sim.TU; %400
bound.TOFGA2_max = 900*3600*24/sim.TU; 
% N REV GA (25)
bound.N_REVGA2_min = 0; %0
bound.N_REVGA2_max = 1; %3
% Incoming Leg
% v_inf_magn GA (26)
bound.v_inf_magnGA2_min = 0;
bound.v_inf_magnGA2_max = 15/sim.DU*sim.TU; % Longusky
% azimuth GAin (27)
bound.azGA2in_min = -pi;
bound.azGA2in_max = pi;
% elevation GAin (28)
bound.elGA2in_min = -pi/2;
bound.elGA2in_max = pi/2;
% Outgoing leg
% azimuth GAout (29)
bound.azGA2out_min = -pi;
bound.azGA2out_max = pi;
% elevation GAout (30)
bound.elGA2out_min = -pi/2;
bound.elGA2out_max = pi/2;

% x = [MJD0,TOF1,TOF2,TOF3,TOF4,NREV,NREV2,NREV3,NREV4,IDP,v_inf_magn,az,el,...
%      CT1,CTa]
bound.lb = [bound.mjd2000_ed, bound.TOF1_min, bound.TOF2_min, ...
    bound.TOFa_min, bound.TOFb_min, bound.N_REV1_min, bound.N_REV2_min, ...
    bound.N_REVa_min, bound.N_REVb_min, bound.IDP_min, bound.IDP2_min, ...
    bound.v_inf_magn_min, bound.az_min, bound.el_min, bound.CT1_min, bound.CTa_min, ...
    bound.TOFGA_min, bound.N_REVGA_min, bound.v_inf_magnGA_min, bound.azGAin_min, ...
    bound.elGAin_min, bound.azGAout_min, bound.elGAout_min, ...
    bound.TOFGA2_min, bound.N_REVGA2_min, bound.v_inf_magnGA2_min, bound.azGA2in_min, ...
    bound.elGA2in_min, bound.azGA2out_min, bound.elGA2out_min]; % Lower bound
bound.ub = [bound.mjd2000_ld, bound.TOF1_max, bound.TOF2_max, ...
    bound.TOFa_max, bound.TOFb_max, bound.N_REV1_max, bound.N_REV2_max, ...
    bound.N_REVa_max, bound.N_REVb_max, bound.IDP_max, bound.IDP2_max, ...
    bound.v_inf_magn_max, bound.az_max, bound.el_max, bound.CT1_max, bound.CTa_max, ...
    bound.TOFGA_max, bound.N_REVGA_max, bound.v_inf_magnGA_max, bound.azGAin_max, ...
    bound.elGAin_max, bound.azGAout_max, bound.elGAout_max, ...
    bound.TOFGA2_max, bound.N_REVGA2_max, bound.v_inf_magnGA2_max, bound.azGA2in_max, ...
    bound.elGA2in_max, bound.azGA2out_max, bound.elGA2out_max]; % Upper bound

%% Constraints
constr.A = []; % linear inequality constraints
constr.b = []; % linear inequality constraints
constr.Aeq = []; % linear equality constraints
constr.beq = []; % linear equality constraints
constr.nonlcon = []; % linear equality constraints
% if you want to restrict x(2) and x(10) to be integers, set IntCon to [2,10].
% ga(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)
% (6) NREV1, (7) NREV2, (8) NREVa, (9) NREVb, (10) IDP, (11) IDP2, (18) NREV GA, (25) NREV GA2
constr.IntCon = [6,7,8,9,10,11,18,25];

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
FitnessFunction = @(x) ff_ea_2SC_LT_RV_soo_NLI(x,sim,data); % Function handle to the fitness function
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

sol.departure_date = mjd20002date(x(1)*sim.TU/(3600*24));
sol.TOF1 = x(2)*sim.TU/86400;
sol.TOF2 = x(3)*sim.TU/86400;
sol.TOFa = x(4)*sim.TU/86400;
sol.TOFb = x(5)*sim.TU/86400;
sol.TOF1_ADIM = x(2);
sol.TOF2_ADIM = x(3);
sol.TOFa_ADIM = x(4);
sol.TOFb_ADIM = x(5);
sol.CT1 = x(15)*sim.TU/86400;
sol.CTa = x(16)*sim.TU/86400;
sol.CT1_ADIM = x(15);
sol.CTa_ADIM = x(16);
sol.mission_duration_d_SC1 = sol.TOF1+sol.CT1+sol.TOF2;
sol.mission_duration_d_SC2 = sol.TOFa+sol.CTa+sol.TOFb;
sol.mission_duration_y_SC1 = sol.mission_duration_d_SC1/365;
sol.mission_duration_y_SC2 = sol.mission_duration_d_SC2/365;
sol.arrival_date_latest = mjd20002date(x(1)*sim.TU/(3600*24)+max(sol.mission_duration_d_SC1,sol.mission_duration_d_SC2));

% sol.mass_fraction = Fval;
sol.Nrev = [x(6), x(7), x(8), x(9)];
sol.v_inf_magn = x(12)*sim.DU/sim.TU;
sol.az = x(13);
sol.el = x(14);

%% characteristic quantities plot
[output, r_encounter, v_encounter, sol] = plot_ff_ea_2SC_LT_GA_RV_soo_NLI(x,sim,data,sol);

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
R_transf_orbit_2 = rotate_local2ecplitic(r_encounter.astD1,r_transf_orbit_2,sim.n_sol,output.Href.leg2);

r_transf_orbit_a  = [output.r.lega.*cos(output.theta.lega), ...
    output.r.lega.*sin(output.theta.lega), output.z.lega];
R_transf_orbit_a = rotate_local2ecplitic(r_encounter.EA,r_transf_orbit_a,sim.n_sol,output.Href.lega);

r_transf_orbit_b  = [output.r.legb.*cos(output.theta.legb), ...
    output.r.legb.*sin(output.theta.legb), output.z.legb];
R_transf_orbit_b = rotate_local2ecplitic(r_encounter.astDa,r_transf_orbit_b,sim.n_sol,output.Href.legb);

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
    '*','Color',colors(8,:),'DisplayName','Earth Dep')
plot3(r_encounter.astA1(1),r_encounter.astA1(2),r_encounter.astA1(3),...
    '^','Color',colors(3,:),'DisplayName',sol.asteroid_1+'Arr')
plot3(r_encounter.astD1(1),r_encounter.astD1(2),r_encounter.astD1(3),...
    '*','Color',colors(3,:),'DisplayName',sol.asteroid_1+'Dep')
plot3(r_encounter.astA2(1),r_encounter.astA2(2),r_encounter.astA2(3),...
    '^','Color',colors(4,:),'DisplayName',sol.asteroid_2+'Arr')
plot3(r_encounter.astAa(1),r_encounter.astAa(2),r_encounter.astAa(3),...
    '^','Color',colors(5,:),'DisplayName',sol.asteroid_a+'Arr')
plot3(r_encounter.astDa(1),r_encounter.astDa(2),r_encounter.astDa(3),...
    '*','Color',colors(5,:),'DisplayName',sol.asteroid_a+'Dep')
plot3(r_encounter.astAb(1),r_encounter.astAb(2),r_encounter.astAb(3),...
    '^','Color',colors(6,:),'DisplayName',sol.asteroid_b+'Arr')
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

%% Checks
output.tGA1 - TOFGA1
output.t1 - TOF1
