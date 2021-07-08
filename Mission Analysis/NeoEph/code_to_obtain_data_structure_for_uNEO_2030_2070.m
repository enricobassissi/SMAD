%% --------------------------------------------------------------------- %%
%% ------- code to finally obtain the stuff for uNEO3 valid ------------ %%
%% ----- in the interval for close approach analysis ------------------- %%
%% -------------------- 2030 - 2070 ------------------------------------ %%
%% --------------------------------------------------------------------- %%
%% ---------- OBTAIN THE VALID ASTEROIDS PART -------------------------- %%
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

%% Python module extraction and elaboration
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(path_py+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

% % Obtain the Sentry risk list
% sentry_risk_names=py.neo_api_function.get_sentry_risk_list();
% % Obtain the NEOCC esa risk list (.txt required, can be download by TODO add link)
% esa_risk_names=py.neo_api_function.extract_esa_name_from_file("esa_risk_list.txt");
% % Merge risk list
% py_risk_list=py.neo_api_function.merge_risk_lists(esa_risk_names, sentry_risk_names);
% 
% % Query JPL SBDB for the bodies in the risk list
% py_dict_risk_list=py.neo_api_function.get_dict(py_risk_list);
% py_ps_vector = module.palermo_scale(py_dict_risk_list);
% ps_vector=cellfun(@double,cell(py_ps_vector));
% 
% % Dictionary extraction of the filtered list
% py_selected_asteroids = py.neo_api_function.refined_selection(py_dict_risk_list);
% 
% % Properties dictionary and name list
% py_selected_asteroids_dict = py_selected_asteroids{1};
% selected_asteroids_names = cellfun(@string,cell(py_selected_asteroids{2}));
% 
% % Magnitude H
% py_MOID_H = py.neo_api_function.MOID_H(py_selected_asteroids_dict);
% H = cellfun(@double,cell(py_MOID_H{2})');
% 
% % Selected Asteroids Characteristics Cell
% [selected_asteroids_orbital_elements_and_sigma, orbital_elements_units] = ...
%     get_orbital_elements_and_sigma(selected_asteroids_names,py_selected_asteroids_dict);
% 
% % initialise
% % a_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (1)
% % e_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (2)
% % incl_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (3)
% % OM_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);  % (4)
% % om_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (5)
% % extraction of data
% idx = 0;
% for i = 1:length(selected_asteroids_orbital_elements_and_sigma)
% %     if selected_asteroids_orbital_elements_and_sigma{i}(1,1) < 2 % pre pruning of semimajor axis
%         idx = idx+1;
%         a_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(1,1);
%         e_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(2,1);
%         incl_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(3,1);
%         OM_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(4,1);
%         om_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(5,1);
%         sel_asteroids_names(idx,1) = selected_asteroids_names(i);
% %     end
% end
% clearvars i idx
% 
% data_elements_matrix = [sel_asteroids_names,a_asteroids,e_asteroids,incl_asteroids,...
%     OM_asteroids,om_asteroids];
% % table of results
% TABLE = table(sel_asteroids_names,H,a_asteroids,e_asteroids,incl_asteroids,OM_asteroids,...
%     om_asteroids);
% 
% %% cut down the asteroids on the orbital parameters
% TF_a_up = str2double(data_elements_matrix(:,2))<1.5; % check a upper bound
% TF_a_low = str2double(data_elements_matrix(:,2))> 0.8; % check a lower bound
% TF_e = str2double(data_elements_matrix(:,3))<0.5; % check e
% TF_i = str2double(data_elements_matrix(:,4))<4; % check i
% 
% data_elements_matrix_cut = data_elements_matrix(and(and(and(TF_a_up,TF_e),TF_a_low),TF_i),:);
% 
% clearvars TF_a_up TF_a_low TF_e TF_i
% 
% %% no pruning permutation
% p_number = 2;
% [PermutationMatrix_noPruning, ~] = permnUnique(sel_asteroids_names, p_number);
% 
% %% Local Pruning on i and w_up
% p_number = 2;
% [asteroid_names, PermutationMatrix_after, HowMany_after] = ...
%             sequences_local_pruning2(data_elements_matrix_cut, p_number);
% 
% %% preparation for the saving to other scripts
% data.asteroid_names = asteroid_names;
% data.PermutationMatrix = PermutationMatrix_after;
% data.HowMany = HowMany_after;
% data.data_elements_matrix = data_elements_matrix_cut;
% data.p_number = p_number;
% 
% %% check on param of the pruned sequence
% for i = 1:HowMany_after
%     for j = 1:p_number
%         idx_ast_considered = find(PermutationMatrix_after(i,j)==data_elements_matrix(:,1));
%         a_perm_mat(i,j) = data_elements_matrix(idx_ast_considered,2);
%         e_perm_mat(i,j) = data_elements_matrix(idx_ast_considered,3);
%         i_perm_mat(i,j) = data_elements_matrix(idx_ast_considered,4);
%     end
% end
% clearvars i j idx_ast_considered
% 
% data.a_perm_mat = a_perm_mat;
% data.e_perm_mat = e_perm_mat;
% data.i_perm_mat = i_perm_mat;

%% --------------------------------------------------------------------- %%
%% -------------- START OF EPHEMERIDES INTERPOLATION PART -------------- %%
%% --------------------------------------------------------------------- %%
%% Asteroid real ephemerides
PointOfView = 'Sun';
epoch_start = '2030-01-01';
epoch_stop = '2070-01-01';
step = '3d';
type_elements = 'elements';

asteroid_names = {'2009TD17';'2011BP40';'2020VV';'2021JE1'};

horizons_data = cell(length(asteroid_names),1);
for name = 1:length(asteroid_names)
    % data extraction section
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(string(asteroid_names(name))),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data{name} = double(py_data); % [a] in AU; [i,OM,om,theta] in deg
    % unwrapping of theta for better interpolation
    % pass to rad all the [i, OM, om theta]
    horizons_data{name}(:,3:6) = deg2rad(horizons_data{name}(:,3:6));
    horizons_data{name}(:,6) = unwrap(horizons_data{name}(:,6));
end
clearvars name

%% Fourier interpolation
t_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),length(horizons_data{1,1}(:,1)));

for name = 1:length(asteroid_names)
    % Fourier interpolation
    N = 400; % number of query point, interpolation points
    V_interp_ft_2{name,1} = interpft(horizons_data{name}(:,1:5),N,1); % 1 by column, 2 by rows
    V_interp_ft_2{name,1}(:,1) = V_interp_ft_2{name,1}(:,1).*AU; % semimajor axis AU -> km
    t_vector2 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N);
    N2 = 1000; % number of query point, interpolation points
    V_interp_ft_22{name,1} = interpft(horizons_data{name}(:,6),N2,1); % 1 by column, 2 by rows
    t_vector22 = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N2);
end
clearvars name

data.y_interp_ft_1_5 = V_interp_ft_2;
data.y_interp_ft_6 = V_interp_ft_22;
data.t_vector_1_5 = t_vector2;
data.t_vector_6 = t_vector22;
data.asteroid_names = asteroid_names;

%% save the data structure
save data_elements_matrix_4_CA_Analysis.mat data
