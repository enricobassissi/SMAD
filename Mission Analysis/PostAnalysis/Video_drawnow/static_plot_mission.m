%% --------------------------------------------------------------------- %%
%% --------------- STATIC PLOT OF THE WHOLE MISSION -------------------- %%
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
str_path=split(pwd, 'PostAnalysis\Video_drawnow');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
dt_path=string(str_path(1))+'TrajOptimisation\direct transcript\follower';
addpath(genpath(dt_path));
fun_path=string(str_path(1))+'TrajOptimisation\direct transcript\functions';
addpath(genpath(fun_path));

%% --------------------------------------------------------------------- %%
%% ------------------------- SC1 --------------------------------------- %%
%% --------------------------------------------------------------------- %%
%% Load of SC1 data
load('data_before_video_sc1.mat')

%% find the new moments where the orbit switch
for i=2:n
    if new_time_vect(i) > MJDA1_dim && new_time_vect(i-1) < MJDA1_dim
        phase_change.astA1 = i+3;
    elseif new_time_vect(i) > MJDD1_dim && new_time_vect(i-1) < MJDD1_dim
        phase_change.astD1 = i-1;
    elseif new_time_vect(i) > MJDA2_dim && new_time_vect(i-1) < MJDA2_dim
        phase_change.astA2 = i-1;
    end
end
    
%%
figure()
% plot3(SC1.uniform.R(:,1),SC1.uniform.R(:,2),SC1.uniform.R(:,3));
% plot3(SC1.uniform.R(1:sim.n_sol,1),SC1.uniform.R(1:sim.n_sol,2),SC1.uniform.R(1:sim.n_sol,3),...
%     'Color',colors(1,:),'DisplayName','Cruise - Leg 1');
axis equal; grid on; hold on;
hl1 = plot3( yplot(1:phase_change.astA1,1), yplot(1:phase_change.astA1,2), yplot(1:phase_change.astA1,3),...
    'Color',colors(1,:),'DisplayName','SC1 - Cruise Leg 1');
hr1 = plot3( yplot(phase_change.astA1:phase_change.astD1,1), yplot(phase_change.astA1:phase_change.astD1,2),...
    yplot(phase_change.astA1:phase_change.astD1,3),'Color',colors(1,:));
hr1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hr1.Color(4) = 0.5;
hl2 = plot3( yplot(phase_change.astD1:phase_change.astA2,1), yplot(phase_change.astD1:phase_change.astA2,2),...
    yplot(phase_change.astD1:phase_change.astA2,3),'Color',colors(1,:));
hl2.Annotation.LegendInformation.IconDisplayStyle = 'off';
hl2.Color(4) = 0.5;
hr2 = plot3( yplot(phase_change.astA2:end,1), yplot(phase_change.astA2:end,2),...
    yplot(phase_change.astA2:end,3),'Color',colors(1,:));
hr2.Annotation.LegendInformation.IconDisplayStyle = 'off';
hr2.Color(4) = 0.5;
% plot3(SC1.uniform.R(sim.n_sol+1:2*sim.n_sol,1),SC1.uniform.R(sim.n_sol+1:2*sim.n_sol,2),...
%     SC1.uniform.R(sim.n_sol+1:2*sim.n_sol,3),...
%     'Color',colors(2,:),'DisplayName','Rendezvous - 2020VV');
% plot3(SC1.uniform.R(2*sim.n_sol+1:3*sim.n_sol,1),SC1.uniform.R(2*sim.n_sol+1:3*sim.n_sol,2),...
%     SC1.uniform.R(2*sim.n_sol+1:3*sim.n_sol,3),...
%     'Color',colors(3,:),'DisplayName','Cruise - Leg 2');
% plot3(SC1.uniform.R(3*sim.n_sol+1:4*sim.n_sol,1),SC1.uniform.R(3*sim.n_sol+1:4*sim.n_sol,2),...
%     SC1.uniform.R(3*sim.n_sol+1:4*sim.n_sol,3),...
%     'Color',colors(4,:),'DisplayName','Rendezvous - 2009TD17');
plot3(new_r_encounter.EA(1), new_r_encounter.EA(2), new_r_encounter.EA(3),'*',...
    'Color',colors(8,:),'DisplayName','Earth Dep');
plot3(new_r_encounter.astA1(1), new_r_encounter.astA1(2), new_r_encounter.astA1(3),'^',...
    'Color',colors(3,:),'DisplayName','2020VV Arr');
plot3(new_r_encounter.astD1(1), new_r_encounter.astD1(2), new_r_encounter.astD1(3),'*',...
    'Color',colors(3,:),'DisplayName','2020VV Dep');
plot3(new_r_encounter.astA2(1), new_r_encounter.astA2(2), new_r_encounter.astA2(3),'^',...
    'Color',colors(4,:),'DisplayName','2009TD17 Arr');
plot3(new_r_encounter.astD2(1), new_r_encounter.astD2(2), new_r_encounter.astD2(3),'*',...
    'Color',colors(4,:),'DisplayName','End Mission');
hsun = plot3(0,0,0,'*','Color',colors(4,:));
hsun.Annotation.LegendInformation.IconDisplayStyle = 'off';

[~, h_earth_whole] = plot_object_orbit(new_time_vect(1),'earth',365,sim,data,100,colors,8);
[~, h_mars_whole] = plot_object_orbit(new_time_vect(1),'mars',687,sim,data,100,colors,6);
h_earth_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
h_mars_whole.Annotation.LegendInformation.IconDisplayStyle = 'off';
h_earth_whole.Color(4) = 0.2;
h_mars_whole.Color(4) = 0.2;
legend('show');

xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
