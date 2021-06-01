%% --------------------------------------------------------------------- %%
%% -------------- PERMUTATIONS AND LOCAL BRANCHING PRUNING ------------- %%
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
selected_asteroids_names = cellfun(@string,cell(py_selected_asteroids{2}));

% Magnitude H
py_MOID_H = py.neo_api_function.MOID_H(py_selected_asteroids_dict);
H = cellfun(@double,cell(py_MOID_H{2})');

% Selected Asteroids Characteristics Cell
[selected_asteroids_orbital_elements_and_sigma, orbital_elements_units] = ...
    get_orbital_elements_and_sigma(selected_asteroids_names,py_selected_asteroids_dict);

% initialise
% a_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (1)
% e_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (2)
% incl_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (3)
% OM_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1);  % (4)
% om_asteroids = zeros(length(selected_asteroids_orbital_elements_and_sigma),1); % (5)
% extraction of data
idx = 0;
for i = 1:length(selected_asteroids_orbital_elements_and_sigma)
%     if selected_asteroids_orbital_elements_and_sigma{i}(1,1) < 2 % pre pruning of semimajor axis
        idx = idx+1;
        a_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(1,1);
        e_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(2,1);
        incl_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(3,1);
        OM_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(4,1);
        om_asteroids(idx,1) = selected_asteroids_orbital_elements_and_sigma{i}(5,1);
        sel_asteroids_names(idx,1) = selected_asteroids_names(i);
%     end
end
clearvars i idx

data_elements_matrix = [sel_asteroids_names,a_asteroids,e_asteroids,incl_asteroids,...
    OM_asteroids,om_asteroids];
% table of results
TABLE = table(sel_asteroids_names,H,a_asteroids,e_asteroids,incl_asteroids,OM_asteroids,...
    om_asteroids)

%% cut down the asteroids on the orbital parameters
TF_a_up = str2double(data_elements_matrix(:,2))<1.8; % check a upper bound
TF_a_low = str2double(data_elements_matrix(:,2))> 0.6; % check a lower bound
TF_e = str2double(data_elements_matrix(:,3))<0.7; % check e
TF_i = str2double(data_elements_matrix(:,4))<5; % check i

data_elements_matrix_cut = data_elements_matrix(and(and(and(TF_a_up,TF_e),TF_a_low),TF_i),:);

clearvars TF_a_up TF_a_low TF_e TF_i

%% no pruning permutation
p_number = 4;
[PermutationMatrix_noPruning, ~] = permnUnique(asteroid_names, p_number);

%% Local Pruning on i and w_up
p_number = 4;
[asteroid_names, PermutationMatrix_after, HowMany_after] = ...
            sequences_local_pruning2(data_elements_matrix_cut, p_number);

%% preparation for the saving to other scripts
data.asteroid_names = asteroid_names;
data.PermutationMatrix = PermutationMatrix_after;
data.HowMany = HowMany_after;
data.data_elements_matrix = data_elements_matrix_cut;
data.p_number = p_number;

%% check on param of the pruned sequence
for i = 1:HowMany_after
    for j = 1:p_number
        idx_ast_considered = find(PermutationMatrix_after(i,j)==data_elements_matrix(:,1));
        a_perm_mat(i,j) = data_elements_matrix(idx_ast_considered,2);
        e_perm_mat(i,j) = data_elements_matrix(idx_ast_considered,3);
        i_perm_mat(i,j) = data_elements_matrix(idx_ast_considered,4);
    end
end
clearvars i j idx_ast_considered

data.a_perm_mat = a_perm_mat;
data.e_perm_mat = e_perm_mat;
data.i_perm_mat = i_perm_mat;