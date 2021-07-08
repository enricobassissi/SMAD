%% --------------------------------------------------------------------- %%
%% ---------------------- Try Impulsive mission ------------------------ %%
%% -------------------- to go to those other objects ------------------- %%
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
          0    174  157;... % (12) PURE TEAL
          0    0    0]./255; % (13) BLACK

%% add path of functions and python stuff
str_path=split(pwd, 'MissionExtension');
addpath(genpath(pwd));
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

load('data_MissionExt_SC1.mat')

%% Mission to the one most promising
idx_ast_best_performing = [1,3,4,7,10,12,15,16]; % ehmmm yes, visual look
names_ast_best_performing = cad_objects_2009TD17(idx_ast_best_performing);
params_ast_best_performing = cad_params_2009TD17(idx_ast_best_performing,:);

%% Adimensionalisation
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0

%% --- TO BE COMPLETED, HEURISTIC METHOD
% for i=1:length(names_ast_best_performing)
%     % --- Build the soo and Run the Optimisation --- %
%     % --- Boundaries
%     % Departure dates (1)
%     bound.mjd2000_ed = pystr2mjd2000(date_of_ca(i))*0.9;
%     bound.mjd2000_ld = pystr2mjd2000(date_of_ca(i))*1.1;
%     % TOF1 (2)
%     bound.TOF1_min = 0; % days
%     bound.TOF1_max = 1*365; % days
% 
%     % x = [MJD0,TOF1]
%     bound.lb = [bound.mjd2000_ed, bound.TOF1_min]; % Lower bound
%     bound.ub = [bound.mjd2000_ld, bound.TOF1_max]; % Upper bound
% 
%     % --- Options
%     options = optimoptions('particleswarm');
%     % options.HybridFcn = @fmincon;
%     options.SwarmSize = 100; % Default is min(100,10*nvars),
%     options.MaxIterations = 30; %  Default is 200*nvars
%     options.MaxStallIterations = 30; % Default 20
%     options.Display = 'iter';
%     options.FunctionTolerance = 1e-6;
% 
%     % % Parallel pool
%     % % Open the parallel pool
%     % par_pool = gcp; 
%     % if isempty(par_pool)
%     %     poolsize = 0;
%     % else
%     %     poolsize = par_pool.NumWorkers;
%     % end
%     % 
%     % options.UseParallel = true;
% 
%     FitnessFunction = @(x) ff_ME_impulsive_soo(x, data, sim); % Function handle to the fitness function
%     numberOfVariables = length(bound.ub); % Number of decision variables
% 
%     tic
%     [x,Fval,exitFlag,Output] = particleswarm(FitnessFunction,numberOfVariables...
%         ,bound.lb,bound.ub,options);
%     el_time_min_pp = toc/60;
% end

%% --- TRIPLE LOOP METHOD, OVER ASTEROID & DEPARTURE & TOF
% --- POSITIONS
% Asteroids
mjd2000_start = pystr2mjd2000(end_of_SC1_mission);
for i=1:length(date_of_ca)
    temp_max_date_ca(i) = pystr2mjd2000(date_of_ca(i));
end
mjd2000_stop = 1.2*max(temp_max_date_ca);
epoch_stop = mjd20002pystr(mjd2000_stop);
step = '2d';
type_elements = 'Vectors';
base_time_vec = [mjd2000_start:2:mjd2000_stop]';

% 2009TD17
py_TL_2009TD17 = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),py.str(PointOfView),...
                  py.str(end_of_SC1_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
TL_2009TD17 = double(py_TL_2009TD17); % [x,y,z] in AU; [vx,vy,vz] in AU/day

% --- data extraction section
% complete orbit
for name=1:length(names_ast_best_performing)
    py_best_ast_data = py.neo_api_function.get_horizons_ephemerides(py.str(names_ast_best_performing(name)),py.str(PointOfView),...
                  py.str(end_of_SC1_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_best_ast_data{name} = double(py_best_ast_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day

%     figure()
%     plot3(horizons_best_ast_data{name}(:,1),horizons_best_ast_data{name}(:,2),horizons_best_ast_data{name}(:,3))
%     legend('names_ast_best_performing(name)')
end


%% TRIPLE LOOP
N = 200;
% max_TOF = 0.2*365;
max_TOF = 5*(max(temp_max_date_ca) - mjd2000_start);
time_to_go_vect = linspace(0, max_TOF, 2*N); % days
departure_mjd2000_vec = linspace(mjd2000_start,1.01*max(temp_max_date_ca),N); % mjd2000
tic
for k=1:length(names_ast_best_performing) % asteroid arrival selected
    ast_to_go = names_ast_best_performing(k);
    TL_ast = horizons_best_ast_data{k};
    for i=1:length(departure_mjd2000_vec) % departure variation
        day_dep = departure_mjd2000_vec(i);
        states_ast_1 = interp1(base_time_vec,TL_2009TD17,day_dep)'; % [x,y,z] in AU; [vx,vy,vz] in AU/day
        states_ast_1(4:6) = states_ast_1(4:6)/86400*sim.TU; % it's adim now
        for j=1:length(time_to_go_vect) % tof variation
            tof = time_to_go_vect(j);
            states_ast_2 = interp1(base_time_vec,TL_ast,day_dep+tof)';
            states_ast_2(4:6) = states_ast_2(4:6)/86400*sim.TU; % it's adim now
            [dv,~,~]=lambert_solver_rendezvous(states_ast_1(1:3),states_ast_2(1:3),...
                states_ast_1(4:6),states_ast_2(4:6),tof*86400/sim.TU,sim.mu); % everything adim... AU, adim vel, adim mu
            dvtot(i,j,k) = dv*sim.DU/sim.TU; % ri adimensionalise -> km/s
%             states_ast_2(4:6) = states_ast_2(4:6)/86400*sim.DU;
%             [dvtot(i,j,k),~,~]=lambert_solver_rendezvous(states_ast_1(1:3)*sim.DU,states_ast_2(1:3)*sim.DU,...
%                 states_ast_1(4:6),states_ast_2(4:6),tof*86400,sim.mu_dim);
        end
    end
end
el_time = toc;
clearvars i j k
for i=1:length(departure_mjd2000_vec)
    time_vec_contour(i) = datenum(datetime(mjd20002date(departure_mjd2000_vec(i))));
end
clearvars i
colors2(:,:,1) = [30   51   120;...% (2) TRUSTY AZURE +2
          0    97   158;...% (3) TRUSTY AZURE +1
          0    155  219;...
          0    142  122;...% (6) PURE TEAL +1
          0    174  157]./255; % (7) PURE TEAL
colors2(:,:,2) = [255  204  78; ...%
          251  171  24;... % (7) ENLIGHT YELLOW
          244  121  32;... % (8) ENLIGHT YELLOW +1
          207  29   57;... % (9) EXCITE RED +1
          150  1    54]./255;    % (10) EXCITE RED +2
% figure()
% for i=1:2
% %     figure()
%     contour(time_vec_contour,time_to_go_vect,dvtot(:,:,i)',[0:0.1:3]);%,[0:0.1:1],'ShowText','on'
%     AX = gca;
%     cololors2 = colors2(:,:,i);
%     colormap(AX,cololors2)
%     hcb = colorbar;
%     hcb.Title.String = "\Delta V [km/s]";
%     xlabel('Date Dep'); ylabel('TOF [d]');
%     ax = gca;
%     ax.XTick=time_vec_contour(1:10:end) ;
%     xtickangle(30)
%     datetick('x','yyyy mmm dd','keepticks')
%     xlim ([time_vec_contour(1) time_vec_contour(end)])
%     hold on   
% end

% --- plot the first 2 that are the best performing overall
figure()
tiledlayout(2,1)
ax1 = nexttile;
contour(time_vec_contour,time_to_go_vect,dvtot(:,:,1)',[0:0.1:2]);
cololors2 = colors2(:,:,1);
colormap(ax1,cololors2)
hcb = colorbar;
hcb.Title.String = "\Delta V [km/s]";
xlabel('Date Dep'); ylabel('TOF [d]');
ax = gca;
ax.XTick=time_vec_contour(1:10:end) ;
xtickangle(30)
datetick('x','yyyy mmm dd','keepticks')
xlim ([time_vec_contour(1) time_vec_contour(end)])
title(names_ast_best_performing(1))
xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca(1))))),'--','LineWidth',2)

ax2 = nexttile;
contour(time_vec_contour,time_to_go_vect,dvtot(:,:,2)',[0:0.1:3]);
cololors2 = colors2(:,:,2);
colormap(ax2,cololors2)
hcb = colorbar;
hcb.Title.String = '\Delta V [km/s]';
xlabel('Date Dep'); ylabel('TOF [d]');
ax = gca;
ax.XTick=time_vec_contour(1:10:end);
xtickangle(30)
datetick('x','yyyy mmm dd','keepticks')
xlim ([time_vec_contour(1) time_vec_contour(end)])
title(names_ast_best_performing(2))
xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca(2))))),'--','LineWidth',2)

