%% MOVIE PLOT OF THE WHOLE MISSION - old version
%{
and the full orbit of inner solar system in that same 
time of the mission
%}
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
% str_path_1=split(pwd, 'follower');
% imp_path=string(str_path_1(1))+'functions';
% addpath(genpath(imp_path));

%% load and textures
load('Mission_160kg_dry.mat')
load('MA5_160kg_dry.mat')
load('data_elements_matrix_44_63_2SC.mat')

%% the time vector
% MA.SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
%     SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000
time_vect = MA.SC1.mjd2000_dep+MA.SC1.uniform.time/86400;

time_vect 

% -- length
n = length(time_vect);

%% sc trajectory
% MA.SC1.uniform.R = [SC1.leg1.EPS.R_cartesian; SC1.coasting.leg1.r_ast; 
%     SC1.leg2.EPS.R_cartesian; SC1.coasting.leg2.r_ast;]./sim.DU; % AU
yplot = MA.SC1.uniform.R/sim.DU;

% % % -- coasting between the two legs
% [SC1.coasting.leg1] = relative_position_coasting_stuff(SC1.asteroid_1,SC1.mjd2000_dep,...
%     SC1.TOF1,SC1.CT1,data.n_int,EphData,sim,colors);%leg1.timead(end)*sim.TU/86400


%% full orbit propagations
% --- same time of propagations between planets/asteroids and sc
for k=1:n
    % -- mercury
	[kepME,ksun] = uplanet(time_vect(k), 1);
	[rME(k,:), vME(k,:)] = sv_from_coe(kepME,ksun);  
% 	V1(k,:)=vME;
    % -- venus
	[kepVE,ksun] = uplanet(time_vect(k), 2);
	[rVE(k,:), vVE] = sv_from_coe(kepVE,ksun);  
% 	V2(k,:)=vVE;
    % -- Earth
	[kepE,ksun] = uplanet(time_vect(k), 3);
	[rEA(k,:), vEA] = sv_from_coe(kepE,ksun);  
% 	V3(k,:)=vE;
    % -- mars
	[kepMA,ksun] = uplanet(time_vect(k), 4);
	[rMA(k,:), vMA] = sv_from_coe(kepMA,ksun);  
% 	V4(k,:)=vMA;
    
    % --- ast 1
    [kep_ast1] = uNEO3(time_vect(k),MA.SC1.asteroid_1,data);
    [r_ast1(k,1:3),v_ast1] = sv_from_coe(kep_ast1,ksun); % km, km/s
    % --- ast 2
    [kep_ast2] = uNEO3(time_vect(k),MA.SC1.asteroid_2,data);
    [r_ast2(k,1:3),v_ast2] = sv_from_coe(kep_ast2,ksun); % km, km/s
end

% -- adimensionalise
rME=rME/sim.DU; rVE=rVE/sim.DU; rEA=rEA/sim.DU; rMA=rMA/sim.DU;
r_ast1 = r_ast1/sim.DU; r_ast2 = r_ast2/sim.DU;

%% sc trajectory propagation
% --- limits of the camera pov
lim_max_helio = max(vecnorm(rMA,2,2));
rel_sc_ast1 = yplot-r_ast1;
% lim_max_ast1 = max(vecnorm(rel_sc_ast1,2,2));
lim_max_ast1 =0.15;
min_rel_sc_ast1_non_zero = 10;
for i = 1:n
    if norm(rel_sc_ast1(i,:)) > 0
        temp_nz = norm(rel_sc_ast1(i,:));
        if temp_nz < min_rel_sc_ast1_non_zero
            min_rel_sc_ast1_non_zero = temp_nz;
        end
    end
end
rel_sc_ast2 = yplot-r_ast2;
min_rel_sc_ast2_non_zero = 10;
lim_max_ast2 = max(vecnorm(rel_sc_ast2,2,2));
for i = 1:n
    if norm(rel_sc_ast2(i,:)) > 0
        temp_nz = norm(rel_sc_ast2(i,:));
        if temp_nz < min_rel_sc_ast2_non_zero
            min_rel_sc_ast2_non_zero = temp_nz;
        end
    end
end

% --- movie
figure('Name','Mission Video')
% figure('renderer','openGL','position',[2    32   800   663]);
set(gca,'nextplot','replacechildren');
v = VideoWriter('mission1','MPEG-4');
% v.Quality = 95;
% axis([-lim_max  lim_max     -lim_max   lim_max     -lim_max   lim_max ])
open(v);  
axis equal
k=0;
% set(gca, 'CameraPosition', [0 lim_max lim_max]);
for i=1:n
    if i <= 95 || and(i > 205, i <= 295)
        % --- normal condition, nominal trajectory heliocentric
        h0 = plot3( yplot(i,1), yplot(i,2), yplot(i,3),'*-','Color','g','MarkerSize',2);
        hold on  
%         full_orbit_repr('mercury')
        h1= plot3( rME(i,1), rME(i,2), rME(i,3),'*-','Color','k','MarkerSize',2);
%         full_orbit_repr('venus')
        h2= plot3( rVE(i,1), rVE(i,2), rVE(i,3),'*-','Color','y','MarkerSize',2);
%         full_orbit_repr('earth')
        h3= plot3( rEA(i,1), rEA(i,2), rEA(i,3),'*-','Color','b','MarkerSize',2);
%         full_orbit_repr('mars')
        h4= plot3( rMA(i,1), rMA(i,2), rMA(i,3),'*-','Color','r','MarkerSize',2);

        view([lim_max_helio lim_max_helio lim_max_helio])  
        k=k+1;
    elseif i>95 && i<=205%norm(rel_sc_ast1(i,1:3)) <= min_rel_sc_ast1_non_zero
        % --- we are in rendezvous with first ast -> pov of the ast
        h01 = plot3( rel_sc_ast1(i,1), rel_sc_ast1(i,2), rel_sc_ast1(i,3),'*-','Color','g','MarkerSize',2);
        hold on  
%         full_orbit_repr('mercury')
%         h1= plot3( R1(i,1)./AU, R1(i,2)./AU, R1(i,3)./AU,'*-','Color','k','MarkerSize',2);
%         full_orbit_repr('venus')
%         h2= plot3( R2(i,1)./AU, R2(i,2)./AU, R2(i,3)./AU,'*-','Color','y','MarkerSize',2);
%         full_orbit_repr('earth')
%         h3= plot3( R3(i,1)./AU, R3(i,2)./AU, R3(i,3)./AU,'*-','Color','b','MarkerSize',2);
%         full_orbit_repr('mars')
%         h4= plot3( R4(i,1)./AU, R4(i,2)./AU, R4(i,3)./AU,'*-','Color','r','MarkerSize',2);
%         full_orbit_repr('jupiter')
%         h5= plot3( R5(i,1)./AU, R5(i,2)./AU, R5(i,3)./AU,'*-','Color',[255,140,0]./255,'MarkerSize',2);

        view([lim_max_ast1 lim_max_ast1 lim_max_ast1])  
        k=k+1;
    elseif i>295 %norm(rel_sc_ast2(i,1:3)) <= min_rel_sc_ast2_non_zero
        % --- we are in rendezvous with second ast -> pov of the ast
        h02 = plot3( rel_sc_ast2(i,1), rel_sc_ast2(i,2), rel_sc_ast2(i,3),'*-','Color','g','MarkerSize',2);
%         full_orbit_repr('saturn')
%         h6= plot3( R6(k,1)./AU, R6(k,2)./AU, R6(k,3)./AU,'*-','Color',[217, 83, 25]./255,'MarkerSize',2);
%         full_orbit_repr('uranus')
%         h7= plot3( R7(k,1)./AU, R7(k,2)./AU, R7(k,3)./AU,'*-','Color','c','MarkerSize',2);
%         full_orbit_repr('neptune')
%         h8= plot3( R8(k,1)./AU, R8(k,2)./AU, R8(k,3)./AU,'*-','Color','b','MarkerSize',2);    

        view([lim_max_ast2 lim_max_ast2 lim_max_ast2])
        k=k+1;
    
    end

    drawnow

    zlabel('AU')    
    grid minor
%     date_view=mjd20002date(ode_data_plot(i)/(3600*24));
%     title (sprintf('date of evaluation %d %d %d', date_view(1:3)));
%     set(get(gca,'title'),'Position',[15 13 1.00011])

    frame = getframe(gcf);
    writeVideo(v,frame);
	if i <= 95 || and(i > 205, i <= 295)
        delete(h0);
        delete(h1);
        delete(h2);
        delete(h3);
        delete(h4);
    elseif i>95 && i<=205
        delete(h0);
        delete(h01);
	elseif i>295
        delete(h0);
        delete(h02);
	end

end

close(v);
