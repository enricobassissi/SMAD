%% --------------------------------------------------------------------- %%
%% ---------------------- close approach analysis ---------------------- %%
%% ----------------- for payload detection from ground ----------------- %%
%% --------------------------------------------------------------------- %%
%% Default options
clear; close all; clc;

set(0, 'DefaultTextFontSize', 20) % modify it if too small
set(0, 'DefaultAxesFontSize', 20) % modify it if too small
set(0, 'DefaultLegendFontSize', 20) % modify it if too small
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.8)

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

%% work environment setup
path_str=split(pwd, 'PostAnalysis\SimoneDellAgnello');
path_util=string(path_str(1))+'Utils';
addpath(genpath(path_util));
path_py=string(path_str(1))+'PyInterface\NEO_API_py';
addpath(genpath(path_py));
% dt_path=string(str_path(1))+'TrajOptimisation\direct transcript\follower';
% addpath(genpath(dt_path));
% fun_path=string(str_path(1))+'TrajOptimisation\direct transcript\functions';
% addpath(genpath(fun_path));

load('data_elements_matrix_4_CA_Analysis.mat')

LD = 384399; % km
AU = 149597870.7; % km

%% - earth eph from jpl changes something?
%% Python module extraction and elaboration
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(path_py+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

% Earth 1 year ephemerides
epoch_start = '2030-01-01';
epoch_stop = '2070-01-01';
step = '3d';
type_elements = 'Vectors';

py_data_Earth = py.neo_api_function.get_earth_ephemerides(py.str(epoch_start),...
                py.str(epoch_stop),py.str(step),py.str(type_elements));
data_Earth = double(py_data_Earth);

t_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),length(data_Earth));

%% test to see if it works
% --- closest approach moment under analysis
% 2033-05-14.55420
mjd_ca_2009td17(1) = date2mjd2000([2033,05,14,13,18,0]);

% -- planet and asteroid positions
% Earth
[kep_EA,ksun] = uplanet(mjd_ca_2009td17(1), 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun);

% Asteroid
[kep_ast_Ast] = uNEO3(mjd_ca_2009td17(1),'2009TD17',data); % [km,-,rad,rad,rad,wrapped rad]
[rAst, vAst] = sv_from_coe(kep_ast_Ast,ksun); % km, km/s

vel_rel(1) = sqrt((vAst(1)-vEA(1))^2 + (vAst(2)-vEA(2))^2 + (vAst(3)-vEA(3))^2); % km/s

% -- proj in ast-earth line
% versor of asteroid position, positive outward
r_Ast_EA = rAst - rEA;
r_Ast_EA_vers = r_Ast_EA/norm(r_Ast_EA);
% versor of velocity , clockwise
v_Ast_EA = vAst - vEA;
v_Ast_EA_vers = v_Ast_EA/norm(v_Ast_EA);
% versor normal to the plane of the ast orbit
h_Ast_EA_vers = cross(r_Ast_EA_vers,v_Ast_EA_vers);
% tg versor, normal to r_vers and h_vers
tg_Ast_EA_vers = cross(r_Ast_EA_vers,h_Ast_EA_vers);

% projection of rel velocity of the ast wrt the earth
dot_prod_r = dot(v_Ast_EA,r_Ast_EA_vers);
dot_prod_th = dot(v_Ast_EA,tg_Ast_EA_vers);
dot_prod_z = dot(v_Ast_EA,h_Ast_EA_vers);
r_Ast_EA = dot_prod_r*r_Ast_EA_vers;
th_Ast_EA = dot_prod_th*tg_Ast_EA_vers;
z_Ast_EA = dot_prod_z*h_Ast_EA_vers;

velocity_Ast_EA = r_Ast_EA+th_Ast_EA+z_Ast_EA; % in the frame Ast_EA

% the first position should be the direction radial, in the conjunction ast-ea
velocity_Ast_EA(1);

%% all the ast and all the close approaches
% --- aSteroids under analysis
asteroids = {'2009TD17';'2011BP40';'2020VV';'2021JE1'};

% --- closest approach moment under analysis
% 2009TD17
date_ca{1,1} = [2033,05,14,13,18,3;
              2039,05,07,13,44,33;
              2052,09,07,2,26,16;
              2058,05,09,20,17,40;
              2064,04,26,11,17,37];

% 2011BP40
date_ca{2,1} = [2035,09,08,22,0,0;%.91886
                2036,03,27,19,0,0;% .74615
                2045,06,25,20,0,0; %.79779
                2046,01,17,0,55,0;%.07950
                2054,06,23,8,0,0;%.29951
                2055,01,22,6,0,0];%.18660

% 2020VV            
date_ca{3,1} = [2033,10,17,10,0,0;%.37289
                2050,07,20,6,0,0;%.25047
                2050,12,31,12,0,0;%.47584
                2056,06,02,5,0,0;%.16696
                2056,10,31,0,0,0];%.06460

% 2021JE1
date_ca{4,1} = [2049,05,15,5,30,0;%.21074
                2050,05,22,20,20,0;%.80446
                2053,10,21,19,30,0;%.70847
                2054,12,10,6,0,0;%.47387
                2055,12,21,4,0,0;%.35995
                2061,04,13,23,50,0;%.99634
                2061,05,31,14,0,0;%.61165
                2062,05,04,3,0,0;%.14181
                2063,04,30,23,0,0];%.99297

% hr = 24*.47057
% min = 60*.2937
% sec = 60*0.6220

mjd_ca = cell(size(date_ca,1),1);
for i=1:size(date_ca,1)
    for j=1:size(date_ca{i},1)
        mjd_ca{i}(j,1) = date2mjd2000(date_ca{i}(j,:));
    end
end
clearvars i j

for i=1:size(date_ca,1) % rows
    for j=1:size(date_ca{i},1) % cols

    % -- planet and asteroid positions
    % Earth
%     [kep_EA,ksun] = uplanet(mjd_ca{i}(j,1), 3);
%     [rEA, vEA] = sv_from_coe(kep_EA,ksun);
    % [x,y,z] in AU; [vx,vy,vz] in AU/day
    earth_state_vector = interp1(t_vector, data_Earth, mjd_ca{i}(j,1), 'spline'); 
    rEA = earth_state_vector(1:3)'*AU; % km
    vEA = earth_state_vector(4:6)'*AU/86400; % km/s

    % Asteroid
    [kep_ast_Ast] = uNEO3(mjd_ca{i}(j,1),asteroids(i),data); % [km,-,rad,rad,rad,wrapped rad]
    [rAst, vAst] = sv_from_coe(kep_ast_Ast,ksun); % km, km/s

    vel_rel_magn{i,1}(j,1) = sqrt((vAst(1)-vEA(1))^2 + (vAst(2)-vEA(2))^2 + (vAst(3)-vEA(3))^2); % km/s

    % -- proj in ast-earth line
    % versor of asteroid position, positive outward
    r_Ast_EA = rAst - rEA;
    r_Ast_EA_vers = r_Ast_EA/norm(r_Ast_EA);
    % versor of velocity , clockwise
    v_Ast_EA = vAst - vEA;
    v_Ast_EA_vers = v_Ast_EA/norm(v_Ast_EA);
    % versor normal to the plane of the ast orbit
    h_Ast_EA_vers = cross(r_Ast_EA_vers,v_Ast_EA_vers);
    % tg versor, normal to r_vers and h_vers
    tg_Ast_EA_vers = cross(r_Ast_EA_vers,h_Ast_EA_vers);

    % projection of rel velocity of the ast wrt the earth
    dot_prod_r = dot(v_Ast_EA,r_Ast_EA_vers);
    dot_prod_th = dot(v_Ast_EA,tg_Ast_EA_vers);
    dot_prod_z = dot(v_Ast_EA,h_Ast_EA_vers);
    r_Ast_EA = dot_prod_r*r_Ast_EA_vers;
    th_Ast_EA = dot_prod_th*tg_Ast_EA_vers;
    z_Ast_EA = dot_prod_z*h_Ast_EA_vers;

    velocity_Ast_EA = r_Ast_EA + th_Ast_EA + z_Ast_EA; % in the polar frame of Ast_EA

    % the first position should be the direction radial, in the conjunction ast-ea
    % the second is the transversal velocity, the one we want
    transv_vel{i,1}(j,1) = velocity_Ast_EA(2);

    end
end


