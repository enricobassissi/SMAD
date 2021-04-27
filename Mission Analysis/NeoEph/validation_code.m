%% --------------- START UP --------------------------------- %%
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
addpath time
addpath Functions

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
% Import the neo_api py library
module = py.importlib.import_module('neo_api_function');

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
selected_asteroids_names = cellfun(@string,cell(py_selected_asteroids{2}));

% Selected Asteroids Characteristics Cell
[selected_asteroids_orbital_elements_and_sigma, orbital_elements_units] = ...
    get_orbital_elements_and_sigma(selected_asteroids_names,py_selected_asteroids_dict);

incl_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);
e_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);
a_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);
for i = 1:length(selected_asteroids_orbital_elements_and_sigma)
    incl_asteroids(i) = selected_asteroids_orbital_elements_and_sigma{i}(3,1);
    e_asteroids(i) = selected_asteroids_orbital_elements_and_sigma{i}(2,1);
    a_asteroids(i) = selected_asteroids_orbital_elements_and_sigma{i}(1,1);
end
clearvars i

TABLE = table(selected_asteroids_names',a_asteroids,e_asteroids,incl_asteroids)

%% ------------------- Analysis on Keplerian Elements --------------------- %%
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

%% ------------------- EXAMPLE TEST -------------------------- %%
%% INTERPOLATION
t_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),length(horizons_data{1,1}(:,1)));

%% Fourier interpolation
N = 200;
y_interp_ft = interpft(horizons_data{1,1},N,1); % 1 by column, 2 by rows
t_vector2 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N);

%% polynomial interpolation
N = 6;
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
legend('real','fourier 200','poly 6')
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

%% ------------------ VALIDATION ------------------------------- %%
%% ------ COMPARISON WITH REAL EPH AND ERROR COMPARISON WITH UPLANET ------- %%
%% Loop with all asteroids and take the worst
t_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),length(horizons_data{1,1}(:,1)));

for name = 1:length(selected_asteroids_names)
    
    % Fourier interpolation
    N = 200; % number of query point, interpolation points
    V_interp_ft{name,1} = interpft(horizons_data{name},N,1); % 1 by column, 2 by rows
    t_vector2 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N);
    
    % interp1 spline of the FT approx, spline as it will be done in uNEO2()
    V_interp_ft_lin{name,1} = interp1(t_vector2,V_interp_ft{name},t_vector, 'spline');

    % error calculation
    V_err{name,1} = horizons_data{name} - V_interp_ft_lin{name};

    % find the max, mean, std error among all the elements, for all the long time span
    V_max_ft200_err{name,1} = max(V_err{name});
    V_mean_ft200_err{name,1} = mean(V_err{name});
    V_std_ft200_err{name,1} = std(V_err{name});
    
end

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

%% order of magnitude of uplanet earth wrt real horizons earth
earth_horizons_2022_2048 = read_web_form_of_JPL_Horizons('earth_horizons_2020_2050.txt');
earth_uplanet = zeros(length(t_vector),6);
for i = 1:length(t_vector)
    earth_uplanet(i,:) = uplanet(t_vector(i),3);
end
% km, rad -> au, deg
earth_uplanet(:,1) = earth_uplanet(:,1)/AU;
earth_uplanet(:,3:6) = rad2deg(earth_uplanet(:,3:6));

err_uplanet = earth_horizons_2022_2048 - earth_uplanet;

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