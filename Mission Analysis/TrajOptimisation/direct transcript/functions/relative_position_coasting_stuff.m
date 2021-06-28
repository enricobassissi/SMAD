function [coasting] = relative_position_coasting_stuff(ast_name,mjd2000_start,TOF,CT,n,EphData,sim,colors)

arr_mjd2000 = mjd2000_start + TOF;
dep_mjd2000 = mjd2000_start + TOF + CT;
coasting = time_and_distance_at_asteroid(ast_name,arr_mjd2000,dep_mjd2000,EphData,sim,n);

% 3D plot of the relative position of ast during coasting
figure('Name','Coasting rel pos ast earth')
plot3(coasting.rel_pos_ast_earth(:,1),coasting.rel_pos_ast_earth(:,2),coasting.rel_pos_ast_earth(:,3),...
    'DisplayName','Ast 1','Color',colors(1,:));
axis equal; hold on;
plot3(0,0,0,'o','DisplayName','Earth','Color',colors(8,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Coasting - 3D rel pos of Ast wrt EA');

% 3D plot of the relative position of earth wrt all ast, during coasting
figure('Name','Coasting rel pos ast earth')
plot3(coasting.rel_pos_earth_ast(:,1),coasting.rel_pos_earth_ast(:,2),coasting.rel_pos_earth_ast(:,3),...
    'DisplayName','Ast 1','Color',colors(1,:));
axis equal; hold on;
plot3(0,0,0,'o','DisplayName','Ast 1','Color',colors(2,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Coasting - 3D rel pos of EA wrt Ast');

% plot of magnitude of the relative distance ast - earth
figure('Name','magnitude of rel pos')
plot(coasting.time-mjd2000_start,coasting.norm_rel_pos_ast_earth,'DisplayName','Ast 1','Color',colors(1,:))
legend('show')
xlabel('time from departure [d]'); ylabel('norm(REL DIST ast wrt EA) [AU]'); 

% actual plot of asteroid and earth, heliocentric
figure('Name','ast ea')
plot3(coasting.r_ast(:,1)./sim.DU,coasting.r_ast(:,2)./sim.DU,coasting.r_ast(:,3)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','ast');
axis equal; hold on;
plot3(coasting.rEA(:,1)./sim.DU,coasting.rEA(:,2)./sim.DU,coasting.rEA(:,3)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(2,:),'DisplayName','EA');
plot3(0,0,0,'o','Color',colors(4,:))
legend('show')
view(2)

% orbital parameters plot of 1st ast to check
figure('Name','orbital params')
subplot(3,2,1)
plot(coasting.time, coasting.coe_ast(:,1),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','a');
legend('show')
subplot(3,2,2)
plot(coasting.time, coasting.coe_ast(:,2),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','e');
legend('show')
subplot(3,2,3)
plot(coasting.time, coasting.coe_ast(:,3),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','i');
legend('show')
subplot(3,2,4)
plot(coasting.time, coasting.coe_ast(:,4),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','OM');
legend('show')
subplot(3,2,5)
plot(coasting.time, coasting.coe_ast(:,5),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','om');
legend('show')
subplot(3,2,6)
plot(coasting.time, coasting.coe_ast(:,6),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','theta');
legend('show')

end