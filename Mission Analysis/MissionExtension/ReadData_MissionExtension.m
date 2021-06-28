%% --------------------------------------------------------------------- %%
%% ------------------- Read data from Python --------------------------- %%
%% ----------------- for close approach with other objects ------------- %%
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

%% aaao 
name = '2009TD17';
PointOfView = 'Sun';
epoch_start = '2032-12-01'; % end of mission on SC 2, with CoastingTime on 2nd ast = CT1 1st ast

py_data_2009TD17 = module.get_horizons_ephemerides_elements(py.str(name),py.str(PointOfView),...
    py.str(epoch_start));
data_2009TD17 = double(py_data_2009TD17);
jd0 = data_2009TD17(1,7); % date of the required orbit as before
jdf = jd0 + 100; % must be max 1 year more than jd0
n_object_requested = 100;
distance_within = 0.05; % I suppose in AU, 0.0026 it's approx 1 LD
py_CAD_2009TD17 = module.get_close_approach_to_asteroid(py_data_2009TD17,py.float(jd0),py.float(jdf),py.int(n_object_requested),py.float(distance_within));

%%
N = str2double(string(CAD_2009TD17{'count'}));
cad_objects_2009TD17 = [];
cad_params_2009TD17 =  [];
idx = 0;
for i=0:N-1
    [py_CAD_params] = module.get_CAD_params(CAD_2009TD17,py.int(i));
    obj = string(py_CAD_params{1});
    params = double(py_CAD_params{2});
    if obj ~= '2009 TD17'
        idx = idx + 1;
        cad_objects_2009TD17 = [cad_objects_2009TD17; obj];
        cad_params_2009TD17(idx,1:6) = params;
    end 
end
%{
output of cad_params:
mjd of encounter
distance close approach km
relative velocity km/s
absolute magnitude H
orbit condition code
pha flag true -> 1
%}

%% TODO DADADA
%{
PIGLIA STI OGGETTI, VEDI L'ORBITA, IL PUNTO DI CONTATTO E QUANDO ECC

PROPAGA LA DISTANZA RELATIVA TRA I DUE PER VEDERE COME SI EVOLVE

FAI UNA MISSIONCINA PER VEDERE SE SI PUò ANDARE DA UNO ALL'ALTRO E QUANTO
COSTEREBBE
%}