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

%% Analysis on Keplerian Elements
% Asteroids
PointOfView = 'Sun';
epoch_start = '2021-01-01';
epoch_stop = '2026-01-01';
step = '5d';
type_elements = 'ephemerides';

horizons_data = cell(length(selected_asteroids_names),1);
for name = 1:length(selected_asteroids_names)
    % data extraction section
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(selected_asteroids_names(name)),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data{name} = double(py_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day
    
end