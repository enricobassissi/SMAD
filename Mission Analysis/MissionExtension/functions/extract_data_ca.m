function [date_of_ca,cad_objects_2009TD17,cad_params_2009TD17] = extract_data_ca(module,last_ast_name,end_of_SC1_mission)

%% SC1 - Mission Extension Proposal
% --- retrieve the data about the close approach on the final orbit of the 1st SC
% last_ast_name = '2009TD17';
PointOfView = 'Sun';
% end_of_SC1_mission = '2032-12-01'; % end of mission on SC 2, with CoastingTime on 2nd ast = CT1 1st ast

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

% %% Orbit characterisation
% % % Earth 1 year ephemerides
% % PointOfView = 'Sun';
% % epoch_stop = mjd20002pystr(pystr2mjd2000(end_of_SC1_mission)+365);
% % step = '5d';
% % type_elements = 'Vectors';
% % 
% % py_data_Earth = py.neo_api_function.get_earth_ephemerides(py.str(end_of_SC1_mission),...
% %                 py.str(epoch_stop),py.str(step),py.str(type_elements));
% % data_Earth = double(py_data_Earth);
% 
% % Asteroids
% epoch_stop = mjd20002pystr(pystr2mjd2000(end_of_SC1_mission)+2*365);
% step = '5d';
% type_elements = 'Vectors';
% 
% % 2009TD17
% py_orbit_2009TD17 = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),py.str(PointOfView),...
%                   py.str(end_of_SC1_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
% orbit_2009TD17 = double(py_orbit_2009TD17); % [x,y,z] in AU; [vx,vy,vz] in AU/day
% 
% horizons_data = cell(length(cad_objects_2009TD17),1);
for name = 1:length(cad_objects_2009TD17)
    % --- date of the close approach
    start_date_of_encounter = mjd20002pystr(mjd2mjd2000(cad_params_2009TD17(name,1)));
    date_of_ca(name,1) = start_date_of_encounter;
%     end_date_of_encounter = mjd20002pystr(mjd2mjd2000(cad_params_2009TD17(name,1))+5); % 5 days more, it's the step, so that we get only 2 points
%     % --- nominal last asteoroid 
%     py_data_of_encounter = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),...
%         py.str(PointOfView),py.str(start_date_of_encounter),py.str(end_date_of_encounter),...
%         py.str(step),py.str(type_elements));
%     encounter_data_2009TD17_temp = double(py_data_of_encounter);  
%     encounter_data_2009TD17 = encounter_data_2009TD17_temp(1,:); % only the actual date
%     
%     % --- data extraction section
%     % complete orbit
%     py_data = py.neo_api_function.get_horizons_ephemerides(py.str(cad_objects_2009TD17(name)),py.str(PointOfView),...
%                       py.str(end_of_SC1_mission),py.str(epoch_stop),py.str(step),py.str(type_elements));
%     horizons_data{name} = double(py_data); % [x,y,z] in AU; [vx,vy,vz] in AU/day
%     % --- moment of close approach
%     py_data_of_encounter = py.neo_api_function.get_horizons_ephemerides(py.str(cad_objects_2009TD17(name)),...
%         py.str(PointOfView),py.str(start_date_of_encounter),py.str(end_date_of_encounter),...
%         py.str(step),py.str(type_elements));
%     encounter_data_temp = double(py_data_of_encounter);  
%     encounter_data = encounter_data_temp(1,:); % only the actual date
end

end