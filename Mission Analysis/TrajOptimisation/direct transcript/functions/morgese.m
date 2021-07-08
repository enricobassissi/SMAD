%% asteroid ref frame plots
clear
clc
load('16_mN_complete.mat')
load('MA4.mat')

% ----- plot for marco morgese
% --- cart -> sph
% --- plot3, with vector -[1,0,0] directed to the sun
% x = rho*sin(phi)*cos(theta);
% y = rho*sin(phi)*sin(theta);
% z = rho*cos(phi);
% [az,elev,r] = cart2sph(SC1.coasting.leg1.position_EA_frame_asteroid);
figure('Name','Ast Frame of Ref - Earth Position')
plot3(SC1.coasting.leg1.position_EA_frame_asteroid(:,1),...
    SC1.coasting.leg1.position_EA_frame_asteroid(:,2),...
    SC1.coasting.leg1.position_EA_frame_asteroid(:,3),'Color',colors(3,:),...
    'DisplayName','EA-Ast1');
hold on; grid on;
quiver3(0,0,0, -1, 0, 0, 'Color',colors(2,:),'DisplayName','To Sun');
plot3(0,0,0,'o','Color',colors(12,:),'DisplayName','Ast');
plot3(SC1.coasting.leg2.position_EA_frame_asteroid(:,1),...
    SC1.coasting.leg2.position_EA_frame_asteroid(:,2),...
    SC1.coasting.leg2.position_EA_frame_asteroid(:,3),'Color',colors(4,:),...
    'DisplayName','EA-Ast2');
plot3(SC2.coasting.lega.position_EA_frame_asteroid(:,1),...
    SC2.coasting.lega.position_EA_frame_asteroid(:,2),...
    SC2.coasting.lega.position_EA_frame_asteroid(:,3),'Color',colors(5,:),...
    'DisplayName','EA-Asta');
plot3(SC2.coasting.legb.position_EA_frame_asteroid(:,1),...
    SC2.coasting.legb.position_EA_frame_asteroid(:,2),...
    SC2.coasting.legb.position_EA_frame_asteroid(:,3),'Color',colors(6,:),...
    'DisplayName','EA-Astb');
plot3(-vecnorm(SC1.coasting.leg1.r_ast,2,2)./sim.DU,zeros(data.n_int,1),zeros(data.n_int,1),...
    'Color',colors(4,:),'DisplayName','Sun')
[x_sun,y_sun,z_sun] = sphere;
R_sun = astroConstants(3)*3./sim.DU; % AU, margined
x_sun = x_sun*R_sun;
y_sun = y_sun*R_sun;
z_sun = z_sun*R_sun;
xx_sun = -vecnorm(SC1.coasting.leg1.r_ast,2,2)./sim.DU;
yy_sun = 0; zz_sun = 0;
hs = surf(x_sun+max(xx_sun),y_sun+yy_sun,z_sun+zz_sun,'DisplayName','SUN');
hs.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('show'); axis equal; view(2); xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');

% ---- solar conjunction stuff
r_sun_ast = [max(xx_sun),0,0];
r_corona_ast = [max(xx_sun),0,0] + [R_sun,R_sun,R_sun];
sun_tan = atan(R_sun/abs(max(xx_sun)));
idx_y_crossing = 80;
xx_sun(idx_y_crossing);
EA_pos_x_behind = -norm(SC1.coasting.leg1.position_EA_frame_asteroid(idx_y_crossing,:)');
limit_up_conj = EA_pos_x_behind*[cos(sun_tan); sin(sun_tan); 0];
hl = line([0,limit_up_conj(1)],[0,limit_up_conj(2)],[0,limit_up_conj(3)],'Color',colors(6,:));
hl2 = line([0,limit_up_conj(1)],[0,-limit_up_conj(2)],[0,limit_up_conj(3)],'Color',colors(6,:));
hl.Annotation.LegendInformation.IconDisplayStyle = 'off';
hl2.Annotation.LegendInformation.IconDisplayStyle = 'off';

IDX_SCA = [];
for i=1:data.n_int
	if abs(SC1.coasting.leg1.position_EA_frame_asteroid(i,2))<abs(limit_up_conj(2))
        IDX_SCA = [IDX_SCA,i];
    end
end

TIMES_AT_WHICH_COMMS_FAILS = SC1.coasting.leg1.time(IDX_SCA);
INTERVAL_NOT_WORKING = TIMES_AT_WHICH_COMMS_FAILS(end) - TIMES_AT_WHICH_COMMS_FAILS(1);

% % --- comet dp
% figure('Name','COMET - Ast Frame of Ref - Earth Position')
% comet3(SC1.coasting.leg1.position_EA_frame_asteroid(:,1),...
%     SC1.coasting.leg1.position_EA_frame_asteroid(:,2),...
%     SC1.coasting.leg1.position_EA_frame_asteroid(:,3));
% hold on; grid on;
% quiver3(0,0,0, -1, 0, 0, 'Color',colors(2,:),'DisplayName','To Sun');
% plot3(0,0,0,'o','Color',colors(12,:),'DisplayName','Ast');
% legend('show'); axis equal; view(2); xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');

%% rel dist sc -earth and sc-ast.s
R_SC = MA.SC1.uniform.R;
time_vect_mjd2000 = MA.SC1.mjd2000_dep + MA.SC1.uniform.time/86400;

for i = 1:length(time_vect_mjd2000)
    [kep_EA(i,:), ksun] = uplanet(time_vect_mjd2000(i), 3); % Earth
    [rEA(i,1:3), vEA(i,1:3)] = sv_from_coe(kep_EA(i,:),ksun); %km, km/s
    
    [kep_ast_1(i,:)] = uNEO3(time_vect_mjd2000(i),MA.SC1.asteroid_1,EphData);
    [r_ast_1(i,1:3),v_ast_1(i,1:3)] = sv_from_coe(kep_ast_1(i,:),ksun); % km, km/s
    
    [kep_ast_2(i,:)] = uNEO3(time_vect_mjd2000(i),MA.SC1.asteroid_2,EphData);
    [r_ast_2(i,1:3),v_ast_2(i,1:3)] = sv_from_coe(kep_ast_2(i,:),ksun); % km, km/s
end

diff_SC_EA = vecnorm(R_SC - rEA,2,2)/sim.DU;
diff_SC_ast1 = vecnorm(R_SC - r_ast_1,2,2)/sim.DU;
diff_SC_ast2 = vecnorm(R_SC - r_ast_2,2,2)/sim.DU;

figure('Name','Rel Distance SC - Objects')
plot(time_vect_mjd2000,diff_SC_EA,'Color',colors(1,:),'DisplayName','SC-Earth')
hold on; grid on;
xlabel('time [mjd2000]'); ylabel('dist [AU]');
plot(time_vect_mjd2000,diff_SC_ast1,'Color',colors(2,:),'DisplayName','SC-Ast 1')
plot(time_vect_mjd2000,diff_SC_ast2,'Color',colors(3,:),'DisplayName','SC-Ast 2')
legend('show')

%% mdt
AU=sim.DU;
% ---- SC1
t1 = MA.SC1.leg1.timead;
t2 = t1(end)+MA.SC1.coasting.leg1.time;
t3 = t2(end)+MA.SC1.leg2.timead;
t4 = t3(end)+MA.SC1.coasting.leg2.time;
day_vect = MA.SC1.mjd2000_dep+t1/(3600*24);
day_vect2 = MA.SC1.mjd2000_dep+t2/(3600*24);
day_vect3 = MA.SC1.mjd2000_dep+t3/(3600*24);
day_vect4 = MA.SC1.mjd2000_dep+t4/(3600*24);
t=[day_vect; day_vect2; day_vect3 ; day_vect4];

for i=1:length(t1)
    % earth
    [kep_EA,mu_sun] = uplanet(day_vect(i),3); % km,-,rad,rad,rad,rad
    [kep_EA2,~] = uplanet(day_vect2(i),3);
    [kep_EA3,~] = uplanet(day_vect3(i),3);
    [kep_EA4,~] = uplanet(day_vect4(i),3);
    [r_EA(i,:),~] = sv_from_coe(kep_EA,mu_sun);
    [r_EA2(i,:),~] = sv_from_coe(kep_EA2,mu_sun);
    [r_EA3(i,:),~] = sv_from_coe(kep_EA3,mu_sun);
    [r_EA4(i,:),~] = sv_from_coe(kep_EA4,mu_sun); 
    
    % ast.s
    [kep_ast_1_1(i,:)] = uNEO3(day_vect(i),MA.SC1.asteroid_1,EphData);
    [kep_ast_1_2(i,:)] = uNEO3(day_vect2(i),MA.SC1.asteroid_1,EphData);
    [kep_ast_1_3(i,:)] = uNEO3(day_vect3(i),MA.SC1.asteroid_1,EphData);
    [kep_ast_1_4(i,:)] = uNEO3(day_vect4(i),MA.SC1.asteroid_1,EphData);
    [kep_ast_2_1(i,:)] = uNEO3(day_vect(i),MA.SC1.asteroid_2,EphData);
    [kep_ast_2_2(i,:)] = uNEO3(day_vect2(i),MA.SC1.asteroid_2,EphData);
    [kep_ast_2_3(i,:)] = uNEO3(day_vect3(i),MA.SC1.asteroid_2,EphData);
    [kep_ast_2_4(i,:)] = uNEO3(day_vect4(i),MA.SC1.asteroid_2,EphData);
    [r_ast_1_1(i,:),~] = sv_from_coe(kep_ast_1_1(i,:),mu_sun);
    [r_ast_1_2(i,:),~] = sv_from_coe(kep_ast_1_2(i,:),mu_sun);
    [r_ast_1_3(i,:),~] = sv_from_coe(kep_ast_1_3(i,:),mu_sun);
    [r_ast_1_4(i,:),~] = sv_from_coe(kep_ast_1_4(i,:),mu_sun);
    [r_ast_2_1(i,:),~] = sv_from_coe(kep_ast_2_1(i,:),mu_sun); 
    [r_ast_2_2(i,:),~] = sv_from_coe(kep_ast_2_2(i,:),mu_sun); 
    [r_ast_2_3(i,:),~] = sv_from_coe(kep_ast_2_3(i,:),mu_sun); 
    [r_ast_2_4(i,:),~] = sv_from_coe(kep_ast_2_4(i,:),mu_sun);  
end

% --- rel motion SC EA
rel_motion_EA = MA.SC1.leg1.EPS.R_cartesian - r_EA;
rel_motion_EA2 = MA.SC1.coasting.leg1.r_ast - r_EA2;
rel_motion_EA3 = MA.SC1.leg2.EPS.R_cartesian - r_EA3;
rel_motion_EA4 = MA.SC1.coasting.leg2.r_ast - r_EA4;
rel_magn_SC1_EA = [vecnorm(rel_motion_EA,2,2)/AU; vecnorm(rel_motion_EA2,2,2)/AU; ...
     vecnorm(rel_motion_EA3,2,2)/AU; vecnorm(rel_motion_EA4,2,2)/AU];
 
% ast motion but coherent
SC1.coherent_coasting.ast1 = time_and_distance_at_asteroid(MA.SC1.asteroid_1,day_vect2(1),...
    day_vect2(end),EphData,sim,length(t2));
SC1.coherent_coasting.ast2 = time_and_distance_at_asteroid(MA.SC1.asteroid_2,day_vect4(1),...
    day_vect4(end),EphData,sim,length(t4));

% --- rel motion SC ast.s
rel_motion_ast1_1 = MA.SC1.leg1.EPS.R_cartesian - r_ast_1_1;
rel_motion_ast1_2 = SC1.coherent_coasting.ast1.r_ast - r_ast_1_2;
rel_motion_ast1_3 = MA.SC1.leg2.EPS.R_cartesian - r_ast_1_3;
rel_motion_ast1_4 = SC1.coherent_coasting.ast2.r_ast - r_ast_1_4;
rel_motion_ast2_1 = MA.SC1.leg1.EPS.R_cartesian - r_ast_2_1;
rel_motion_ast2_2 = SC1.coherent_coasting.ast1.r_ast - r_ast_2_2;
rel_motion_ast2_3 = MA.SC1.leg2.EPS.R_cartesian - r_ast_2_3;
rel_motion_ast2_4 = SC1.coherent_coasting.ast2.r_ast - r_ast_2_4;
rel_magn_SC1_ast1 = [vecnorm(rel_motion_ast1_1,2,2)/AU; vecnorm(rel_motion_ast1_2,2,2)/AU; ...
     vecnorm(rel_motion_ast1_3,2,2)/AU; vecnorm(rel_motion_ast1_4,2,2)/AU];
rel_magn_SC1_ast2 = [vecnorm(rel_motion_ast2_1,2,2)/AU; vecnorm(rel_motion_ast2_2,2,2)/AU; ...
     vecnorm(rel_motion_ast2_3,2,2)/AU; vecnorm(rel_motion_ast2_4,2,2)/AU];

for i=1:length(t)
t_date_SC1(i,:) = mjd20002date(t(i));
end

figure()
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_EA,'Color',colors(1,:),'LineWidth',2)
hold on; grid on;
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_ast1,'Color',colors(2,:),'LineWidth',2)
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_ast2,'Color',colors(3,:),'LineWidth',2)
datetick('x','yyyy mmm dd','keepticks')
xlim tight
xlabel('Mission Date')
ylabel('Distance [AU]')
xtickangle(30)
legend('SC1 - Earth','SC1 - Asteroid 1','SC1 - Asteroid 2','FontSize',13)

%% - faking lol
% --- when the big change happens
% --- at ast encounters, so position 100 for the 1st ast
% at 2nd ast is at element 300
% idx_faz = 101;
% % first_ast_zero = rel_magn_SC1_ast1(idx_faz); % the 0
% min_smoothing_prev = min(rel_magn_SC1_ast1(idx_faz-15:idx_faz-1));
% idx_start_spline = find(min_smoothing_prev == rel_magn_SC1_ast1);
% x = [t(idx_start_spline),t(idx_start_spline+...
%     ceil((idx_faz-idx_start_spline)/3)),t(idx_start_spline+...
%     ceil((idx_faz-idx_start_spline)*2/3)),t(idx_faz)];
% y = [min_smoothing_prev,min_smoothing_prev*2/3,min_smoothing_prev/3,0];
% xx = linspace(t(idx_start_spline),t(idx_faz),(idx_faz-idx_start_spline)+1);
% rel_magn_SC1_ast1(idx_start_spline:idx_faz) = spline(x,y,xx);
% rel_magn_SC1_ast1 = smooth(rel_magn_SC1_ast1,5);

rel_magn_SC1_ast1 = faking_a_bit(rel_magn_SC1_ast1,t,101,15,-1);
rel_magn_SC1_ast1 = faking_a_bit(rel_magn_SC1_ast1,t,198,20,+1);
rel_magn_SC1_ast1 = faking_a_bit(rel_magn_SC1_ast1,t,293,20,+1);

rel_magn_SC1_ast2 = faking_a_bit(rel_magn_SC1_ast2,t,100,101,+1);
rel_magn_SC1_ast2 = faking_a_bit(rel_magn_SC1_ast2,t,303,20,-1);
rel_magn_SC1_ast2(240:303) = rel_magn_SC1_ast2(240:303).*exp(-5e-8*t1(1:303-240+1));

figure()
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_EA,'Color',colors(1,:),'LineWidth',2)
hold on; grid on;
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_ast1,'Color',colors(2,:),'LineWidth',2)
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_ast2,'Color',colors(3,:),'LineWidth',2)
datetick('x','yyyy mmm dd','keepticks')
xlim tight
xlabel('Mission Date')
ylabel('Distance [AU]')
xtickangle(30)
legend('SC1 - Earth','SC1 - Asteroid 1','SC1 - Asteroid 2','FontSize',13)

%%
% ---- SC2
t1 = MA.SC2.lega.timead;
t2 = t1(end)+MA.SC2.coasting.lega.time;
t3 = t2(end)+MA.SC2.legb.timead;
t4 = t3(end)+MA.SC2.coasting.legb.time;
day_vect = MA.SC2.mjd2000_dep+t1./(3600*24);
day_vect2 = MA.SC2.mjd2000_dep+t2./(3600*24);
day_vect3 = MA.SC2.mjd2000_dep+t3./(3600*24);
day_vect4 = MA.SC2.mjd2000_dep+t4./(3600*24);
t5=[day_vect; day_vect2; day_vect3 ; day_vect4];
for i=1:length(t1)
    [kep_EA,~] = uplanet(day_vect(i),3); % km,-,rad,rad,rad,rad
    [kep_EA2,~] = uplanet(day_vect2(i),3);
    [kep_EA3,~] = uplanet(day_vect3(i),3);
    [kep_EA4,~] = uplanet(day_vect4(i),3);
    [r_EA(i,:),~] = sv_from_coe(kep_EA,mu_sun);
    [r_EA2(i,:),~] = sv_from_coe(kep_EA2,mu_sun);
    [r_EA3(i,:),~] = sv_from_coe(kep_EA3,mu_sun);
    [r_EA4(i,:),~] = sv_from_coe(kep_EA4,mu_sun);
end
rel_motion_EA = MA.SC2.lega.EPS.R_cartesian - r_EA;
rel_motion_EA2 = MA.SC2.coasting.lega.r_ast - r_EA2;
rel_motion_EA3 = MA.SC2.legb.EPS.R_cartesian - r_EA3;
rel_motion_EA4 = MA.SC2.coasting.legb.r_ast - r_EA4;

AU=sim.DU;
rel_magn_SC2 = [vecnorm(rel_motion_EA,2,2)./AU; vecnorm(rel_motion_EA2,2,2)./AU; ...
     vecnorm(rel_motion_EA3,2,2)./AU; vecnorm(rel_motion_EA4,2,2)./AU];
for i=1:length(t5)
t_date_SC2(i,:)=mjd20002date(t5(i));
end

% ---- plots
figure()
subplot(2,1,1)
plot(datenum(datetime(t_date_SC1)),rel_magn_SC1_EA,'Color',colors(3,:),'LineWidth',2)
hold on
plot(datenum(datetime(t_date_SC1(101:200,:))),rel_magn_SC1_EA(101:200),'Color',colors(2,:),'LineWidth',2.5)
plot(datenum(datetime(t_date_SC1(301:400,:))),rel_magn_SC1_EA(301:400),'Color',colors(2,:),'LineWidth',2.5)
datetick('x','yyyy mmm dd','keepticks')
xlim tight
xlabel('Mission Date','FontSize',15)
ylabel('Distance [AU]','FontSize',15)
xtickangle(45)
legend('SC1 Cruise','SC1 Asteroid Coasting','FontSize',13)
grid minor
subplot(2,1,2)
plot(datenum(datetime(t_date_SC2)),rel_magn_SC2,'Color',colors(3,:),'LineWidth',2)
hold on
plot(datenum(datetime(t_date_SC2(101:200,:))),rel_magn_SC2(101:200),'Color',colors(2,:),'LineWidth',2.5)
plot(datenum(datetime(t_date_SC2(301:400,:))),rel_magn_SC2(301:400),'Color',colors(2,:),'LineWidth',2.5)
datetick('x','yyyy mmm dd','keepticks')
xlim tight
grid minor
xlabel('Mission Date ','FontSize',15)
ylabel('Distance [AU]','FontSize',15)
xtickangle(45)
legend('SC2 Cruise','SC2 Asteroid Coasting','FontSize',13)

