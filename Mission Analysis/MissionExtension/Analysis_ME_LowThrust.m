%% Low Thrust Mission Extension
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
me_path=string(str_path(1))+'MissionExtension';
addpath(genpath(me_path));
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

%% Adimensionalisation
sim.mu_dim    = 132712440018              ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7               ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0

%% Append multiple 100 days pieces of cloase approach
% structure:
% 1 jd of encounter
% 2 distance close approach km
% 3 relative velocity km/s
% 4 absolute magnitude H
% 5 orbit condition code
% 6 pha flag true -> 1

% --- segment 1
last_ast_name = '2009TD17';
segment_1_date = '2032-12-01';
[date_of_ca_1,cad_objects_2009TD17_1,cad_params_2009TD17_1] = extract_data_ca(module,last_ast_name,...
    segment_1_date);
% --- segment 2
segment_2_date = mjd20002pystr(pystr2mjd2000(segment_1_date)+100);
[date_of_ca_2,cad_objects_2009TD17_2,cad_params_2009TD17_2] = extract_data_ca(module,last_ast_name,...
    segment_2_date);
% --- segment 3 
segment_3_date = mjd20002pystr(pystr2mjd2000(segment_2_date)+100);
[date_of_ca_3,cad_objects_2009TD17_3,cad_params_2009TD17_3] = extract_data_ca(module,last_ast_name,...
    segment_3_date);
% --- segment 4 --- and we got more than 1 year
segment_4_date = mjd20002pystr(pystr2mjd2000(segment_3_date)+100);
[date_of_ca_4,cad_objects_2009TD17_4,cad_params_2009TD17_4] = extract_data_ca(module,last_ast_name,...
    segment_4_date);

% --- append togheter
date_of_ca = [date_of_ca_1;date_of_ca_2;date_of_ca_3;date_of_ca_4];
cad_objects_2009TD17 = [cad_objects_2009TD17_1;cad_objects_2009TD17_2;...
    cad_objects_2009TD17_3;cad_objects_2009TD17_4];
cad_params_2009TD17 = [cad_params_2009TD17_1;cad_params_2009TD17_2;...
    cad_params_2009TD17_3;cad_params_2009TD17_4];

%% Lexell
% data
Lexell_id = find(cad_objects_2009TD17 == 'Lexell');
Lexell_cad = cad_params_2009TD17(Lexell_id,:);
Lexell_date_encounter = mjd2date(Lexell_cad(1));

%% to cut it off
TF = cad_objects_2009TD17 == 'Lexell';
cad_objects_2009TD17 = cad_objects_2009TD17(~TF);
date_of_ca = date_of_ca(~TF);
cad_params_2009TD17 = cad_params_2009TD17(~TF,:);

%% call to Analysis of the Mission Extension LT
end_of_SC1_mission = segment_1_date;
N = 100;
sim.M_end = 160; % kg
sim.Isp = 3200/sim.TU; %s -> adim
sim.TOF_imposed_flag = 1;
sim.direction = 1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.n_sol = N;                      % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 
 
[dvtot_temp,time_vec_contour,time_to_go_vect,base_time_vec,...
    horizons_best_ast_data_temp,TL_2009TD17] = triple_loop_lowthrust_ME(sim,N,...
    end_of_SC1_mission,date_of_ca,last_ast_name,cad_objects_2009TD17);

%% plot the ones with possibility to be done
% --- plot the first 2 that are the best performing overall
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

idx = 0;
for i=1:size(dvtot_temp,3)
    if min(min(dvtot_temp(:,:,i))) < 100 && min(min(dvtot_temp(:,:,i))) ~= 0
        idx = idx+1;
        dvtot_SC1(:,:,idx) = dvtot_temp(:,:,i);
        cad_objects_SC1(idx,1) = cad_objects_2009TD17(i,1);
        date_of_ca_SC1(idx,1) = date_of_ca(i,1);
        horizons_best_ast_data_SC1{idx} = horizons_best_ast_data_temp{i}; % [x,y,z] in AU; [vx,vy,vz] in AU/day
    end
end
clearvars i idx
        
        
for i = 1:2:size(dvtot_SC1,3)
    figure()
    tiledlayout(2,1)
    ax1 = nexttile;
    contour(time_vec_contour,time_to_go_vect,dvtot_SC1(:,:,i)',[0:10:100]);
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
    title(cad_objects_SC1(i))
    
    xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca_SC1(i))))),'--','LineWidth',2)

    ax2 = nexttile;
    contour(time_vec_contour,time_to_go_vect,dvtot_SC1(:,:,i+1)',[0:0.1:2]);
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
    title(cad_objects_SC1(i+1))
    
    xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca_SC1(i+1))))),'--','LineWidth',2)
end
clearvars i