%% --------------------------------------------------------------------- %%
%% -------- VALIDATION CODE FOR EPH MODEL USED FOR ASTEROIDS ----------- %%
%% --------------- COMPARISON BETWEEN REAL EPH FROM JPL AND ------------ %%
%% ----------- INTERPOLATED WITH FOURIER SERIES DATA ------------------- %%
%% --------------------------------------------------------------------- %%
%% Setup for default options
set(0, 'DefaultTextFontSize', 22)
set(0, 'DefaultAxesFontSize', 22)
set(0, 'DefaultLegendFontSize', 22)
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 2)
format short
%% Initialise the environment
clear; close all; clc;
path_str=split(pwd, 'NeoEph');
path_utils=string(path_str(1))+'Utils';
addpath(genpath(path_utils));
path_py=string(path_str(1))+'PyInterface\NEO_API_py';
addpath(genpath(path_py));
path_neoeph=string(path_str(1))+'NeoEph';
addpath(genpath(path_neoeph));

AU = astroConstants(2);
muSun = astroConstants(4);

% Palette ESA
colors = [0    50   71;... % DEEP SPACE
          207  29   57;... % EXCITE RED +1
          0    142  122;... % PURE TEAL +1
          251  171  24;... % ENLIGHT YELLOW
          244  121  32;... % ENLIGHT YELLOW +1
          150  1    54;... % EXCITE RED +2
          167  85   52;... % ENLIGHT YELLOW +2
          0    97   158;... % TRUSTY AZURE +1
          30   51   120;... % TRUSTY AZURE +2
          0    103  98;... % PURE TEAL +2
          51   94   111;... % DEEP SPACE -1
          0    0    0]./255; % BLACK

%% Python module extraction and elaboration
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(path_py+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

% Obtain the Sentry risk list
sentry_risk_names=py.neo_api_function.get_sentry_risk_list();
% Obtain the NEOCC esa risk list (.txt required, can be download by TODO add link)
esa_risk_names=py.neo_api_function.extract_esa_name_from_file("esa_risk_list.txt");
% Merge risk list
py_risk_list=py.neo_api_function.merge_risk_lists(esa_risk_names, sentry_risk_names);

% Query JPL SBDB for the bodies in the risk list
py_dict_risk_list=py.neo_api_function.get_dict(py_risk_list);
py_ps_vector = module.palermo_scale(py_dict_risk_list);
ps_vector=cellfun(@double,cell(py_ps_vector));

% Dictionary extraction of the filtered list
py_selected_asteroids = py.neo_api_function.refined_selection(py_dict_risk_list);

% Properties dictionary and name list
py_selected_asteroids_dict = py_selected_asteroids{1};
selected_asteroids_names_not_filtered = cellfun(@string,cell(py_selected_asteroids{2}));

% Selected Asteroids Characteristics Cell
[selected_asteroids_orbital_elements_and_sigma, orbital_elements_units] = ...
    get_orbital_elements_and_sigma(selected_asteroids_names_not_filtered,py_selected_asteroids_dict);

incl_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);
e_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);
a_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);
for i = 1:length(selected_asteroids_orbital_elements_and_sigma)
    incl_asteroids(i) = selected_asteroids_orbital_elements_and_sigma{i}(3,1);
    e_asteroids(i) = selected_asteroids_orbital_elements_and_sigma{i}(2,1);
    a_asteroids(i) = selected_asteroids_orbital_elements_and_sigma{i}(1,1);
end
clearvars i

TABLE = table(selected_asteroids_names_not_filtered',a_asteroids,e_asteroids,incl_asteroids)

%% --------------------------------------------------------------------- %%
%% ---------------- GET HORIZONS DATA OF THE ASTEROIDS ----------------- %%
%% ------------------------ FROM 2020 TO 2050 -------------------------- %%
%% ------------------ Analysis on Keplerian Elements ------------------- %%
%% --------------------------------------------------------------------- %%
% I actually want only the 42 of the local pruning
load('bbbb.mat');
selected_asteroids_names = aaa;
% selected_asteroids_names = ["2012BK14"];

% Asteroids
PointOfView = 'Sun';
epoch_start = '2020-01-01';
epoch_stop = '2050-01-01';
step = '1d';
type_elements = 'elements';

horizons_data = cell(length(selected_asteroids_names),1);
for name = 1:length(selected_asteroids_names)
    % data extraction section
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(selected_asteroids_names(name)),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data{name} = double(py_data); % [a] in AU; [i,OM,om,theta] in deg
    % unwrapping of theta for better interpolation
    horizons_data{name}(:,6) = rad2deg(unwrap(deg2rad(horizons_data{name}(:,6))));
    
end
clearvars name

%% --------------------------------------------------------------------- %%
%% ----------------------- EXAMPLE TEST -------------------------------- %%
%% --------------------------------------------------------------------- %%
%% INTERPOLATION
t_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),length(horizons_data{1,1}(:,1)));

%% Fourier interpolation
N = 500;
y_interp_ft = interpft(horizons_data{1,1},N,1); % 1 by column, 2 by rows
t_vector2 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N);

%% polynomial interpolation
N = 3;
y_polyfit = polyfit(t_vector,horizons_data{1,1}(:,1),N);
y_polyval = polyval(y_polyfit,t_vector);

%% plots examples
% semimajor axis to see if it works
figure()
plot(t_vector,horizons_data{1,1}(:,1))
hold on
plot(t_vector2,y_interp_ft(:,1))
plot(t_vector,y_polyval)
xlabel('mjd2000'); ylabel('a [AU]');
legend('Horizon','Fourier N200','Poly 3')
title('semimajor axis for',selected_asteroids_names(1))

% true anomaly
figure()
plot(t_vector,horizons_data{1,1}(:,6))
hold on
plot(t_vector2,y_interp_ft(:,6))
xlabel('mjd2000'); ylabel('$\theta$ [deg]');
legend('real','fourier 200')
title('$\theta$ for',selected_asteroids_names(1))

%% interp1 spline
y_interp_ft_lin = interp1(t_vector2,y_interp_ft,t_vector, 'spline');

%% error calculation
err = horizons_data{1,1} - y_interp_ft_lin;

format long
max_ft200_err = max(err);
mean_ft200_err = mean(err);
std_ft200_err = std(err);

%% --------------------------------------------------------------------- %%
%% ------------------ PARAMETERS VALIDATION ---------------------------- %%
%% ----------- COMPARISON FT INTERPOLATION WRT REAL EPH ---------------- %%
%% ---------------- AND ERROR COMPARISON WRT UPLANET ------------------- %%
%% --------------------------------------------------------------------- %%
%% Loop with all asteroids and take the worst
t_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),length(horizons_data{1,1}(:,1)));

for name = 1:length(selected_asteroids_names)
    
    % Fourier interpolation
    N = 400; % number of query point, interpolation points
    V_interp_ft_2{name,1} = interpft(horizons_data{name}(:,1:5),N,1); % 1 by column, 2 by rows
    t_vector2 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N);
    N2 = 1000; % number of query point, interpolation points
    V_interp_ft_22{name,1} = interpft(horizons_data{name}(:,6),N2,1); % 1 by column, 2 by rows
    t_vector22 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N2);
    
    % interp1 spline of the FT approx, spline as it will be done in uNEO2()
    V_interp_ft_lin{name,1}(:,1:5) = interp1(t_vector2,V_interp_ft_2{name},t_vector, 'spline');
    V_interp_ft_lin{name,1}(:,6) = interp1(t_vector22,V_interp_ft_22{name},t_vector, 'spline');

    % error calculation
    V_err{name,1} = horizons_data{name} - V_interp_ft_lin{name};

    % find the max, mean, std error among all the elements, for all the long time span
    V_max_ft200_err{name,1} = max(V_err{name});
    V_mean_ft200_err{name,1} = mean(V_err{name});
    V_std_ft200_err{name,1} = std(V_err{name});
    
end
clearvars name

format long
% max of each error element among all of the asteroids
% initialize as the first error vector. if we find worse, than replace it
VV_max_ft200_err = V_max_ft200_err{1}; 
VV_mean_ft200_err = V_mean_ft200_err{1}; 
VV_std_ft200_err = V_std_ft200_err{1}; 

for i = 1:length(V_max_ft200_err{1}) % orb params
    for j=2:length(V_max_ft200_err) % asteroids cells
        if V_max_ft200_err{j}(i) > V_max_ft200_err{j-1}(i)
            VV_max_ft200_err(i) = V_max_ft200_err{j}(i);
        end
        if V_mean_ft200_err{j}(i) > V_mean_ft200_err{j-1}(i)
            VV_mean_ft200_err(i) = V_mean_ft200_err{j}(i);
        end
        if V_std_ft200_err{j}(i) > V_std_ft200_err{j-1}(i)
            VV_std_ft200_err(i) = V_std_ft200_err{j}(i);
        end
    end
end
clearvars i j

%% order of magnitude of uplanet earth wrt real horizons earth
earth_horizons_2020_2050 = read_web_form_of_JPL_Horizons('earth_horizons_2020_2050.txt');
earth_uplanet = zeros(length(t_vector),6);
for i = 1:length(t_vector)
    earth_uplanet(i,:) = uplanet(t_vector(i),3);
end
clearvars i
% km, rad -> au, deg
earth_uplanet(:,1) = earth_uplanet(:,1)/AU;
earth_uplanet(:,3:6) = rad2deg(earth_uplanet(:,3:6));

err_uplanet = earth_horizons_2020_2050 - earth_uplanet;

format long
max_uplanet_err = max(err_uplanet);
mean_uplanet_err = mean(err_uplanet);
std_uplanet_err = std(err_uplanet);

%% error comparison for validation
err_ratio_max = VV_max_ft200_err/max_uplanet_err
err_ratio_mean = VV_mean_ft200_err/mean_uplanet_err
err_ratio_std = VV_std_ft200_err/std_uplanet_err

fprintf('this fourier transform model interpolation is more precise than the uplanet analytical model \n')
fprintf('if both compared with the result from JPL Horizons Elements output \n')

%% PLOT THE ERROR ON EACH PARAM
% semimajor axis to see if it works
unit_of_measurments = {'a [AU]','e [-]','i [deg]','OM [deg]','om [deg]','th [deg]'};

for i = 1:size(V_err{1},2) % 6 elements
    figure('Name','Error on Elements')
    for j = 1:length(selected_asteroids_names) % 9 asteroids
        h_fig = plot(t_vector',V_err{j}(:,i),'DisplayName',selected_asteroids_names(j));
        hold on
    end
    txt_ylabel = string({'error ', unit_of_measurments(i)});
    xlabel('mjd2000'); ylabel(txt_ylabel);
    legend('show','Location','southeastoutside')
    title('err = real JPL Horizons - FT interpolation')
    saveas(h_fig,sprintf('./Figures/ValidationOfParams/valid_err_param%d.png',i));
end
clearvars i j h_fig

%% error of uplanet plot

idx_th_grater_180 = earth_horizons_2020_2050(:,6)>=180;
for i=1:length(idx_th_grater_180)
    if idx_th_grater_180(i)
        earth_horizons_2020_2050(i,6) = earth_horizons_2020_2050(i,6)-360;
    end
end
    
for i = 1:size(err_uplanet,2)
    figure('Name','Error on Elements of uplanet Earth')
    plot(t_vector',err_uplanet(:,i));
    txt_ylabel = string({'error ', unit_of_measurments(i)});
    xlabel('mjd2000'); ylabel(txt_ylabel);
    legend('show','Location','southeast')
    title('err = real JPL Horizons - uplanet analytical ephemerides')
end
clearvars i

%% true anomaly
% interp1 spline
% V_interp_ft_lin
figure()
plot(t_vector,horizons_data{1,1}(:,6))
hold on
plot(t_vector,V_interp_ft_lin{1,1}(:,6))
xlabel('mjd2000'); ylabel('$\theta$ [deg]');
legend('real','fourier 200')
title('$\theta$ for',selected_asteroids_names(1))

max(horizons_data{1,1}(:,6) - V_interp_ft_lin{1,1}(:,6))
min(horizons_data{1,1}(:,6) - V_interp_ft_lin{1,1}(:,6))
%% --------------------------------------------------------------------- %%
%% ----------------------- ORBIT VALIDATION ---------------------------- %%
%% --------------------------------------------------------------------- %%
% Earth 1 year ephemerides
epoch_start = '2021-01-01';
epoch_stop = '2022-01-01';
step = '5d';
type_elements = 'Vectors';

py_data_Earth = py.neo_api_function.get_earth_ephemerides(py.str(epoch_start),...
                py.str(epoch_stop),py.str(step),py.str(type_elements));
data_Earth = double(py_data_Earth);

% Asteroids
PointOfView = 'Sun';
epoch_start = '2021-01-01';
epoch_stop = '2026-01-01';
step = '5d';
type_elements = 'Vectors';

t_vector_orbit = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),ceil((pystr2mjd2000(epoch_stop)-pystr2mjd2000(epoch_start))/5));

horizons_data_orbit = cell(length(selected_asteroids_names),1);
for name = 1:length(selected_asteroids_names)
    % data extraction section of API data
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(selected_asteroids_names(name)),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data_orbit{name,1} = double(py_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day
    
    % data extraction section of interpolated data
    V_interp_ft_lin_orbit(:,1:5) = interp1(t_vector2,V_interp_ft_2{name},t_vector_orbit, 'spline');
    V_interp_ft_lin_orbit(:,6) = interp1(t_vector22,V_interp_ft_22{name},t_vector_orbit, 'spline');
    % from interp1 they are [a,e,i,OM,om,theta] in [AU,-,deg,deg,deg,deg]
    % but for sv_from_coe I need [a,e,i,OM,om,theta] in [km,-,rad,rad,rad,rad]
    V_interp_ft_lin_orbit(:,1) = V_interp_ft_lin_orbit(:,1)*AU;
    V_interp_ft_lin_orbit(:,3:6) = deg2rad(V_interp_ft_lin_orbit(:,3:6));
    % and re wrapping of theta
    V_interp_ft_lin_orbit(:,6) = wrapping_to_2pi(V_interp_ft_lin_orbit(:,6));
    for i = 1:length(V_interp_ft_lin_orbit)
        [r_ast_params_temp(i,1:3),~] = sv_from_coe(V_interp_ft_lin_orbit(i,:),muSun); % km, km/s
    end
    r_ast_interp_spline_orbit{name,1} = r_ast_params_temp./AU; % km-> AU
    
    % plotting section
    h_fig = figure();
    h_EA = plot3(data_Earth(:,1),data_Earth(:,2),data_Earth(:,3),...
           'LineWidth',2.2,"Color",colors(1,:),'DisplayName','Earth');
    hold on
    h_asteroid = plot3(horizons_data_orbit{name}(:,1),horizons_data_orbit{name}(:,2),horizons_data_orbit{name}(:,3),...
                 'LineWidth',2.2,"Color",colors(2,:),'DisplayName',strcat('API-',selected_asteroids_names(name)));
    
    h_asteroid_interp = plot3(r_ast_interp_spline_orbit{name}(:,1),r_ast_interp_spline_orbit{name}(:,2),r_ast_interp_spline_orbit{name}(:,3),...
         'LineWidth',2.2,"Color",colors(3,:),'DisplayName',strcat('INTERP-',selected_asteroids_names(name)));
    axis equal; grid on;
    xlabel('x [AU]'); ylabel('y [AU]');  zlabel('z [AU]'); 
    legend('show','Location',"southeast")
%     saveas(h_fig,sprintf('./Figures/Orbits/orb%d.png',name)); % questo ?? giusto
    %print(h_fig,sprintf('./Figures/Orbits/orb%d.pdf',name),'-dpdf','-bestfit'); 
    %exportgraphics(gca,sprintf('./Figures/Orbits/orb%d.png',name),'ContentType','image');
end

%% error on the overall orbit for the time plotted
% Loop with all asteroids and take the worst
for name = 1:length(selected_asteroids_names)
    
    % error calculation
    V_orbit_err{name,1} = r_ast_interp_spline_orbit{name} - horizons_data_orbit{name}(:,1:3);

    % find the max, mean, std error among all the elements, for all the long time span
    V_max_orbit_err{name,1} = max(V_orbit_err{name});
    V_mean_orbit_err{name,1} = mean(V_orbit_err{name});
    V_std_orbit_err{name,1} = std(V_orbit_err{name});
    
end

format long
% max of each error element among all of the asteroids
% initialize as the first error vector. if we find worse, than replace it
VV_max_orbit_err = V_max_orbit_err{1}; 
VV_mean_orbit_err = V_mean_orbit_err{1}; 
VV_std_orbit_err = V_std_orbit_err{1}; 

for i = 1:length(V_max_orbit_err{1}) % orb params
    for j=2:length(V_max_orbit_err) % asteroids cells
        if V_max_orbit_err{j}(i) > V_max_orbit_err{j-1}(i)
            VV_max_orbit_err(i) = V_max_orbit_err{j}(i);
        end
        if V_mean_orbit_err{j}(i) > V_mean_orbit_err{j-1}(i)
            VV_mean_orbit_err(i) = V_mean_orbit_err{j}(i);
        end
        if V_std_orbit_err{j}(i) > V_std_orbit_err{j-1}(i)
            VV_std_orbit_err(i) = V_std_orbit_err{j}(i);
        end
    end
end

VV_max_orbit_err %AU
VV_mean_orbit_err %AU
VV_std_orbit_err %AU

%% PLOT EARTH WITH HORIZONS AND WITH UPLANET
coe_earth_hor = zeros(size(earth_horizons_2020_2050,1),size(earth_horizons_2020_2050,2));
coe_earth_hor(:,1) = earth_horizons_2020_2050(:,1)*AU;
coe_earth_hor(:,2) = earth_horizons_2020_2050(:,2);
coe_earth_hor(:,3:6) = deg2rad(earth_horizons_2020_2050(:,3:6));

coe_earth_upl = zeros(size(earth_horizons_2020_2050,1),size(earth_horizons_2020_2050,2));
coe_earth_upl(:,1) = earth_uplanet(:,1)*AU;
coe_earth_upl(:,2) = earth_uplanet(:,2);
coe_earth_upl(:,3:6) = deg2rad(earth_uplanet(:,3:6));
for i=1:size(earth_horizons_2020_2050,1)
    earth_plot_horizon_2020_2050(i,:) = sv_from_coe(coe_earth_hor(i,:),muSun);
    earth_plot_uplanet_2020_2050(i,:) = sv_from_coe(coe_earth_upl(i,:),muSun);
end

h_fig = figure();
h_EA_hor = plot3(earth_plot_horizon_2020_2050(:,1),earth_plot_horizon_2020_2050(:,2),earth_plot_horizon_2020_2050(:,3),...
       'LineWidth',2.2,"Color",colors(1,:),'DisplayName','Earth HOR');
hold on
h_EA_upl = plot3(earth_plot_uplanet_2020_2050(:,1),earth_plot_uplanet_2020_2050(:,2),earth_plot_uplanet_2020_2050(:,3),...
       'LineWidth',2.2,"Color",colors(2,:),'DisplayName','Earth uplanet');
legend('show')

%% delete the python file from this directory
delete('neo_api_function.py');