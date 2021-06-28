%% asteroid ref frame plots
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