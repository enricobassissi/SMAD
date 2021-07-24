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

%% SC1 - Mission Extension Proposal
% --- retrieve the data about the close approach on the final orbit of the 1st SC
last_ast_name = '2009TD17';
PointOfView = 'Sun';
end_of_SC1_mission = '2032-12-01'; % end of mission on SC 2, with CoastingTime on 2nd ast = CT1 1st ast

py_data_2009TD17 = module.get_horizons_ephemerides_elements(py.str(last_ast_name),py.str(PointOfView),...
    py.str(end_of_SC1_mission));
data_2009TD17 = double(py_data_2009TD17);
jd0 = data_2009TD17(1,7); % date of the required orbit as before
jdf = jd0 + 100; % must be max 1 year more than jd0
n_object_requested = 100;
distance_within = 0.05; % I suppose in AU, 0.0026 it's approx 1 LD
py_CAD_2009TD17 = module.get_close_approach_to_asteroid(py_data_2009TD17,py.float(jd0),py.float(jdf),py.int(n_object_requested),py.float(distance_within));

%% extract matlab variables
N = str2double(string(py_CAD_2009TD17{'count'}));
cad_objects_2009TD17 = [];
cad_params_2009TD17 = [];
idx = 0;
for i=0:N-1
    py_CAD_params = module.get_CAD_params(py_CAD_2009TD17,py.int(i));
    obj = strrep(string(py_CAD_params{1}),' ',''); % remove spaces
    params = double(py_CAD_params{2});
    if obj ~= last_ast_name % here with spaces
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

%% Orbit characterisation
% Earth 1 year ephemerides
PointOfView = 'Sun';
epoch_stop = mjd20002pystr(pystr2mjd2000(end_of_SC1_mission)+365);
step = '5d';
type_elements = 'Vectors';

py_data_Earth = py.neo_api_function.get_earth_ephemerides(py.str(end_of_SC1_mission),...
                py.str(epoch_stop),py.str(step),py.str(type_elements));
data_Earth = double(py_data_Earth);

% Asteroids
epoch_stop = mjd20002pystr(pystr2mjd2000(end_of_SC1_mission)+2*365);
step = '5d';
type_elements = 'Vectors';

% 2009TD17
py_orbit_2009TD17 = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),py.str(PointOfView),...
                  py.str(end_of_SC1_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
orbit_2009TD17 = double(py_orbit_2009TD17); % [x,y,z] in AU; [vx,vy,vz] in AU/day

horizons_data = cell(length(cad_objects_2009TD17),1);
for name = 1:length(cad_objects_2009TD17)
    % --- date of the close approach
    start_date_of_encounter = mjd20002pystr(mjd2mjd2000(cad_params_2009TD17(name,1)));
    date_of_ca(name,1) = start_date_of_encounter;
    end_date_of_encounter = mjd20002pystr(mjd2mjd2000(cad_params_2009TD17(name,1))+5); % 5 days more, it's the step, so that we get only 2 points
    % --- nominal last asteoroid 
    py_data_of_encounter = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),...
        py.str(PointOfView),py.str(start_date_of_encounter),py.str(end_date_of_encounter),...
        py.str(step),py.str(type_elements));
    encounter_data_2009TD17_temp = double(py_data_of_encounter);  
    encounter_data_2009TD17 = encounter_data_2009TD17_temp(1,:); % only the actual date
    encounter_data_2009TD17_mat(name,:) = encounter_data_2009TD17;
    
    % --- data extraction section
    % complete orbit
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(cad_objects_2009TD17(name)),py.str(PointOfView),...
                      py.str(end_of_SC1_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data{name} = double(py_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day
    % --- moment of close approach
    py_data_of_encounter = py.neo_api_function.get_horizons_ephemerides(py.str(cad_objects_2009TD17(name)),...
        py.str(PointOfView),py.str(start_date_of_encounter),py.str(end_date_of_encounter),...
        py.str(step),py.str(type_elements));
    encounter_data_temp = double(py_data_of_encounter);  
    encounter_data = encounter_data_temp(1,:); % only the actual date
    encounter_data_mat(name,:) = encounter_data;
    
    % plotting section
    h_fig = figure();
    h_EA = plot3(data_Earth(:,1),data_Earth(:,2),data_Earth(:,3),...
           'LineWidth',2.2,"Color",colors(1,:),'DisplayName','Earth');
    hold on
    h_2009TD17 = plot3(orbit_2009TD17(:,1),orbit_2009TD17(:,2),orbit_2009TD17(:,3),...
           'LineWidth',2.2,"Color",colors(2,:),'DisplayName','2009TD17');
    hp_2009TD17 = plot3(encounter_data_2009TD17(1),encounter_data_2009TD17(2),encounter_data_2009TD17(3),...
        'o',"Color",colors(2,:),'DisplayName','CA 2009TD17');
    h_asteroid = plot3(horizons_data{name}(:,1),horizons_data{name}(:,2),horizons_data{name}(:,3),...
                 'LineWidth',2.2,"Color",colors(3,:),'DisplayName',cad_objects_2009TD17(name));
    hp_asteroid = plot3(encounter_data(1),encounter_data(2),encounter_data(3),...
        '^',"Color",colors(3,:),'DisplayName','CA '+cad_objects_2009TD17(name));
    axis equal; grid on;
    xlabel('x [AU]'); ylabel('y [AU]');  zlabel('z [AU]'); 
    legend('show','Location',"northeast");
%     saveas(h_fig,sprintf('./Figures/SC1/orb%d.png',name));
    %print(h_fig,sprintf('./Figures/Orbits/orb%d.pdf',name),'-dpdf','-bestfit'); 
    %exportgraphics(gca,sprintf('./Figures/Orbits/orb%d.png',name),'ContentType','image');
end

%% data about close appraoch rel vel
cad_2009TD17_2010RM80(1,:) = encounter_data_2009TD17_mat(4,4:6);
cad_2009TD17_2010RM80(2,:) = encounter_data_mat(4,4:6);

cad_2009TD17_2011GJ3(1,:) = encounter_data_2009TD17_mat(5,4:6);
cad_2009TD17_2011GJ3(2,:) = encounter_data_mat(5,4:6);

cad_2009TD17_2010RM80_v_rel = sqrt((cad_2009TD17_2010RM80(1,1)-cad_2009TD17_2010RM80(2,1))^2 + ...
            (cad_2009TD17_2010RM80(1,2)-cad_2009TD17_2010RM80(2,2))^2 + ...
            (cad_2009TD17_2010RM80(1,3)-cad_2009TD17_2010RM80(2,3))^2);
cad_2009TD17_2010RM80_v_rel_kms = cad_2009TD17_2010RM80_v_rel*astroConstants(2)/86400;

cad_2009TD17_2011GJ3_v_rel = sqrt((cad_2009TD17_2011GJ3(1,1)-cad_2009TD17_2011GJ3(2,1))^2 + ...
            (cad_2009TD17_2011GJ3(1,2)-cad_2009TD17_2011GJ3(2,2))^2 + ...
            (cad_2009TD17_2011GJ3(1,3)-cad_2009TD17_2011GJ3(2,3))^2);
cad_2009TD17_2011GJ3_v_rel_kms = cad_2009TD17_2011GJ3_v_rel*astroConstants(2)/86400;
