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
nli_path=string(str_path(1))+'TrajOptimisation\LowThrust';
addpath(genpath(nli_path));

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

% from 2 mesi after l'arrivo a 1 anno max
N = 100;
% --- segment 1
last_ast_name = '2009TD17';
segment_1_date = '2033-04-25';
[date_of_ca_1,cad_objects_2009TD17_1,cad_params_2009TD17_1] = extract_data_ca(module,last_ast_name,...
    segment_1_date);
% --- segment 2
segment_2_date = mjd20002pystr(pystr2mjd2000(segment_1_date)+N);
[date_of_ca_2,cad_objects_2009TD17_2,cad_params_2009TD17_2] = extract_data_ca(module,last_ast_name,...
    segment_2_date);
% --- segment 3 
segment_3_date = mjd20002pystr(pystr2mjd2000(segment_2_date)+N);
[date_of_ca_3,cad_objects_2009TD17_3,cad_params_2009TD17_3] = extract_data_ca(module,last_ast_name,...
    segment_3_date);

% --- append togheter
date_of_ca = [date_of_ca_1;date_of_ca_2;date_of_ca_3];
cad_objects_2009TD17 = [cad_objects_2009TD17_1;cad_objects_2009TD17_2;...
    cad_objects_2009TD17_3];
cad_params_2009TD17 = [cad_params_2009TD17_1;cad_params_2009TD17_2;...
    cad_params_2009TD17_3];

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

%% --------------------------------------------------------------------- %%
%% -------------- call to Analysis of the Mission Extension ------------ %%
%% --------------------------- Impulsive ------------------------------- %%
%% --------------------------------------------------------------------- %%
load('data_ME_append_300d_SC1.mat')

end_of_SC1_mission = segment_1_date;

[dvtot_temp,time_vec_contour,time_to_go_vect,base_time_vec,...
    horizons_best_ast_data_temp,TL_2009TD17] = triple_loop_impulsive_ME2(sim,N,...
    end_of_SC1_mission,last_ast_name,cad_objects_2009TD17);

%% plot the ones with possibility to be done
load('data_MissionExt_SC1_new1.mat')

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
    if min(min(dvtot_temp(:,:,i))) < 3
        idx = idx+1;
        dvtot_SC1(:,:,idx) = dvtot_temp(:,:,i);
        cad_objects_SC1(idx,1) = cad_objects_2009TD17(i,1);
        date_of_ca_SC1(idx,1) = date_of_ca(i,1);
        horizons_best_ast_data_SC1{idx} = horizons_best_ast_data_temp{i}; % [x,y,z] in AU; [vx,vy,vz] in AU/day
    end
end
clearvars i idx
        
for i = 1:2:2 % size(dvtot_SC1,3)
    [row,col] = find(min(min(dvtot_SC1(:,:,1))) == dvtot_SC1(:,:,1));
    figure()
    tiledlayout(2,1)
    ax1 = nexttile;
    contour(time_vec_contour,time_to_go_vect,dvtot_SC1(:,:,i)',[0:0.1:3]);
    hold on
    plot(time_vec_contour(row),time_to_go_vect(col),'o','color',colors(6,:),...
        'markersize',10,'displayname','min')
    xl=xline(time_vec_contour(row),'color',colors(6,:));
    yl=yline(time_to_go_vect(col),'color',colors(6,:));
    xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
    yl.Annotation.LegendInformation.IconDisplayStyle = 'off';
    cololors2 = colors2(:,:,1);
    colormap(ax1,cololors2)
    hcb = colorbar;
    hcb.Title.String = "\Delta V [km/s]";
    xlabel('Date Dep'); ylabel('TOF [d]');
    ax = gca;
    ax.XTick=time_vec_contour(1:4:end) ;
    xtickangle(30)
    datetick('x','yyyy mmm dd','keepticks')
%     xlim ([time_vec_contour(1) time_vec_contour(end)])
    xlim ([time_vec_contour(120) time_vec_contour(end)])
    ylim ([time_to_go_vect(160) time_to_go_vect(240)])
    title(cad_objects_SC1(i))
    
%     xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca_SC1(i))))),'--','LineWidth',2)
    
    [row2,col2] = find(min(min(dvtot_SC1(:,:,2))) == dvtot_SC1(:,:,2));
    ax2 = nexttile;
    contour(time_vec_contour,time_to_go_vect,dvtot_SC1(:,:,i+1)',[0:0.1:3]);
    cololors2 = colors2(:,:,2);
    hold on
    plot(time_vec_contour(row2),time_to_go_vect(col2),'o','color',colors(3,:),...
        'markersize',10,'displayname','min')
    xl=xline(time_vec_contour(row2),'color',colors(3,:));
    yl=yline(time_to_go_vect(col2),'color',colors(3,:));
    xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
    yl.Annotation.LegendInformation.IconDisplayStyle = 'off';
    colormap(ax2,cololors2)
    hcb = colorbar;
    hcb.Title.String = '\Delta V [km/s]';
    xlabel('Date Dep'); ylabel('TOF [d]');
    ax = gca;
    ax.XTick=time_vec_contour(1:5:end);
    xtickangle(30)
    datetick('x','yyyy mmm dd','keepticks')
%     xlim ([time_vec_contour(1) time_vec_contour(end)])
    xlim ([time_vec_contour(50) time_vec_contour(120)])
    ylim ([time_to_go_vect(110) time_to_go_vect(240)])
    title(cad_objects_SC1(i+1))
    
%     xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca_SC1(i+1))))),'--','LineWidth',2)
end
clearvars i

% --- how much mass
dv_min_SC1_IMP = min(min(dvtot_SC1(:,:,2)))*1000;
m_wet_SC1_IMP = 160/exp(-dv_min_SC1_IMP/(9.81*270));

% %% --------------------------------------------------------------------- %%
% %% ------------ real rel vel
% %% Orbit characterisation
% % Earth 1 year ephemerides
% PointOfView = 'Sun';
% epoch_stop = mjd20002pystr(pystr2mjd2000(end_of_SC1_mission)+365);
% step = '5d';
% type_elements = 'Vectors';
% 
% py_data_Earth = py.neo_api_function.get_earth_ephemerides(py.str(end_of_SC1_mission),...
%                 py.str(epoch_stop),py.str(step),py.str(type_elements));
% data_Earth = double(py_data_Earth);
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
% for name = 1:length(cad_objects_2009TD17)
%     % --- date of the close approach
%     start_date_of_encounter = mjd20002pystr(mjd2mjd2000(cad_params_2009TD17(name,1)));
%     date_of_ca(name,1) = start_date_of_encounter;
%     end_date_of_encounter = mjd20002pystr(mjd2mjd2000(cad_params_2009TD17(name,1))+5); % 5 days more, it's the step, so that we get only 2 points
%     % --- nominal last asteoroid 
%     py_data_of_encounter = py.neo_api_function.get_horizons_ephemerides(py.str(last_ast_name),...
%         py.str(PointOfView),py.str(start_date_of_encounter),py.str(end_date_of_encounter),...
%         py.str(step),py.str(type_elements));
%     encounter_data_2009TD17_temp = double(py_data_of_encounter);  
%     encounter_data_2009TD17 = encounter_data_2009TD17_temp(1,:); % only the actual date
%     encounter_data_2009TD17_mat(name,:) = encounter_data_2009TD17;
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
%     encounter_data_mat(name,:) = encounter_data;
%     
%     % plotting section
%     h_fig = figure();
%     h_EA = plot3(data_Earth(:,1),data_Earth(:,2),data_Earth(:,3),...
%            'LineWidth',2.2,"Color",colors(1,:),'DisplayName','Earth');
%     hold on
%     h_2009TD17 = plot3(orbit_2009TD17(:,1),orbit_2009TD17(:,2),orbit_2009TD17(:,3),...
%            'LineWidth',2.2,"Color",colors(2,:),'DisplayName','2009TD17');
%     hp_2009TD17 = plot3(encounter_data_2009TD17(1),encounter_data_2009TD17(2),encounter_data_2009TD17(3),...
%         'o',"Color",colors(2,:),'DisplayName','CA 2009TD17');
%     h_asteroid = plot3(horizons_data{name}(:,1),horizons_data{name}(:,2),horizons_data{name}(:,3),...
%                  'LineWidth',2.2,"Color",colors(3,:),'DisplayName',cad_objects_2009TD17(name));
%     hp_asteroid = plot3(encounter_data(1),encounter_data(2),encounter_data(3),...
%         '^',"Color",colors(3,:),'DisplayName','CA '+cad_objects_2009TD17(name));
%     axis equal; grid on;
%     xlabel('x [AU]'); ylabel('y [AU]');  zlabel('z [AU]'); 
%     legend('show','Location',"northeast");
% %     saveas(h_fig,sprintf('./Figures/SC1/orb%d.png',name));
%     %print(h_fig,sprintf('./Figures/Orbits/orb%d.pdf',name),'-dpdf','-bestfit'); 
%     %exportgraphics(gca,sprintf('./Figures/Orbits/orb%d.png',name),'ContentType','image');
% end
% 
% %% data about close appraoch rel vel
% cad_2009TD17_2010RM80(1,:) = encounter_data_2009TD17_mat(4,4:6);
% cad_2009TD17_2010RM80(2,:) = encounter_data_mat(4,4:6);
% 
% cad_2009TD17_2011GJ3(1,:) = encounter_data_2009TD17_mat(5,4:6);
% cad_2009TD17_2011GJ3(2,:) = encounter_data_mat(5,4:6);
% 
% cad_2009TD17_2010RM80_v_rel = sqrt((cad_2009TD17_2010RM80(1,1)-cad_2009TD17_2010RM80(2,1))^2 + ...
%             (cad_2009TD17_2010RM80(1,2)-cad_2009TD17_2010RM80(2,2))^2 + ...
%             (cad_2009TD17_2010RM80(1,3)-cad_2009TD17_2010RM80(2,3))^2);
% cad_2009TD17_2010RM80_v_rel_kms = cad_2009TD17_2010RM80_v_rel*astroConstants(2)/86400;
% 
% cad_2009TD17_2011GJ3_v_rel = sqrt((cad_2009TD17_2011GJ3(1,1)-cad_2009TD17_2011GJ3(2,1))^2 + ...
%             (cad_2009TD17_2011GJ3(1,2)-cad_2009TD17_2011GJ3(2,2))^2 + ...
%             (cad_2009TD17_2011GJ3(1,3)-cad_2009TD17_2011GJ3(2,3))^2);
% cad_2009TD17_2011GJ3_v_rel_kms = cad_2009TD17_2011GJ3_v_rel*astroConstants(2)/86400;

%% call to Analysis of the Mission Extension LT
end_of_SC1_mission = segment_1_date;
sim.M_end = 160; % kg
sim.Isp = 3200/sim.TU; %s -> adim
sim.TOF_imposed_flag = 1;
sim.direction = 1;                     % direction of integration (1 FW, -1 BW), 
                                       % 1 is like imposing wet mass at beginning
sim.n_sol = N;                      % number of computational nodes
sim.x = linspace(0,1,sim.n_sol)';   % 
asteroid_EOL_chosen = ["2010JK1","2010RM80","2020TW"];

[dvtot_temp_LT,time_vec_contour_LT,time_to_go_vect_LT,base_time_vec_LT,...
    horizons_best_ast_data_temp_LT,TL_2009TD17_LT] = triple_loop_lowthrust_ME2(sim,N,...
    end_of_SC1_mission,last_ast_name,asteroid_EOL_chosen);

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
for i=1:size(dvtot_temp_LT,3)
    idx = idx+1;
%     if min(min(real(dvtot_temp_LT(:,:,i)))) < 500 && min(min(real(dvtot_temp_LT(:,:,i)))) ~= 0
    dvtot_SC1_LT(:,:,idx) = real(dvtot_temp_LT(:,:,i));
    cad_objects_SC1_LT(idx,1) = asteroid_EOL_chosen(i);
%     date_of_ca_SC1_LT(idx,1) = date_of_ca(i,1);
    horizons_best_ast_data_SC1_LT{idx} = horizons_best_ast_data_temp_LT{i}; % [x,y,z] in AU; [vx,vy,vz] in AU/day
%     end
    for j=1:size(dvtot_temp_LT,2)
            for k=1:size(dvtot_temp_LT,1)
                if dvtot_SC1_LT(k,j,i) == 0
                    dvtot_SC1_LT(k,j,i) = NaN;
                end
            end
    end
  
end
clearvars i idx j k
        
        
for i = 1:2:size(dvtot_SC1_LT,3)
    figure()
    tiledlayout(2,1)
    ax1 = nexttile;
    contour(time_vec_contour_LT,time_to_go_vect_LT,dvtot_SC1_LT(:,:,i)',[0:0.3:9]);
    cololors2 = colors2(:,:,1);
    colormap(ax1,cololors2)
    hcb = colorbar;
    hcb.Title.String = "\Delta V [km/s]";
    xlabel('Date Dep'); ylabel('TOF [d]');
    ax = gca;
    ax.XTick=time_vec_contour_LT(1:10:end) ;
    xtickangle(30)
    datetick('x','yyyy mmm dd','keepticks')
    xlim ([time_vec_contour_LT(1) time_vec_contour_LT(end)])
    title(cad_objects_SC1_LT(i))
    
%     xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca_SC1_LT(i))))),'--','LineWidth',2)

    ax2 = nexttile;
    contour(time_vec_contour_LT,time_to_go_vect_LT,dvtot_SC1_LT(:,:,i+1)',[0:0.1:5]);
    cololors2 = colors2(:,:,2);
    colormap(ax2,cololors2)
    hcb = colorbar;
    hcb.Title.String = '\Delta V [km/s]';
    xlabel('Date Dep'); ylabel('TOF [d]');
    ax = gca;
    ax.XTick=time_vec_contour_LT(1:10:end);
    xtickangle(30)
    datetick('x','yyyy mmm dd','keepticks')
    xlim ([time_vec_contour_LT(1) time_vec_contour_LT(end)])
    title(cad_objects_SC1_LT(i+1))
    
%     xline(datenum(datetime(mjd20002date(pystr2mjd2000(date_of_ca_SC1_LT(i+1))))),'--','LineWidth',2)
end
clearvars i

%% only 2010RM80
[row,col] = find(min(min(dvtot_SC1_LT(:,:,2))) == dvtot_SC1_LT(:,:,2));
ax2 = figure();
contour(time_vec_contour_LT,time_to_go_vect_LT,dvtot_SC1_LT(:,:,2)',[0:0.1:5],...
    'DisplayName','Porkchop Plot');
hold on
plot(time_vec_contour_LT(row),time_to_go_vect_LT(col),'o','color',colors(3,:),...
    'markersize',10,'displayname','min')
xl=xline(time_vec_contour_LT(row),'color',colors(3,:));
yl=yline(time_to_go_vect_LT(col),'color',colors(3,:));
xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
yl.Annotation.LegendInformation.IconDisplayStyle = 'off';
cololors2 = colors2(:,:,2);
colormap(ax2,cololors2)
hcb = colorbar;
hcb.Title.String = '\Delta V [km/s]';
xlabel('Date Dep'); ylabel('TOF [d]');
ax = gca;
ax.XTick=time_vec_contour_LT(1:2:end);
xtickangle(30)
datetick('x','yyyy mmm dd','keepticks')
xlim ([time_vec_contour_LT(90) time_vec_contour_LT(110)])
ylim([200 250])
title(cad_objects_SC1_LT(2))
legend('show')

% --- how much mass
dv_min_SC1_LT = min(min(dvtot_SC1_LT(:,:,2)))*1000;
m_wet_SC1_LT = 160/exp(-dv_min_SC1_LT/(9.81*3200));