coasting.ast1.arr_mjd2000 = sol.departure_mjd2000 + sol.TOF1;
coasting.ast1.dep_mjd2000 = sol.departure_mjd2000 + sol.TOF1 + sol.CT1; %
coasting.ast1 = time_and_distance_at_asteroid(sol.asteroid_1,coasting.ast1.arr_mjd2000,coasting.ast1.dep_mjd2000,data,sim);

coasting.ast2.arr_mjd2000 = sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + sol.TOF2;
coasting.ast2.dep_mjd2000 = sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + sol.TOF2 + 50; % arbitrary 50 day of prox ops on last ast
coasting.ast2 = time_and_distance_at_asteroid(sol.asteroid_2,coasting.ast2.arr_mjd2000,coasting.ast2.dep_mjd2000,data,sim);

coasting.asta.arr_mjd2000 = sol.departure_mjd2000 + sol.TOFa;
coasting.asta.dep_mjd2000 = sol.departure_mjd2000 + sol.TOFa + sol.CTa;
coasting.asta = time_and_distance_at_asteroid(sol.asteroid_a,coasting.asta.arr_mjd2000,coasting.asta.dep_mjd2000,data,sim);

coasting.astb.arr_mjd2000 = sol.departure_mjd2000 + sol.TOFa + sol.CTa + sol.TOFb;
coasting.astb.dep_mjd2000 = sol.departure_mjd2000 + sol.TOFa + sol.CTa + sol.TOFb + 50; % arbitrary 50 day of prox ops on last ast
coasting.astb = time_and_distance_at_asteroid(sol.asteroid_b,coasting.astb.arr_mjd2000,coasting.astb.dep_mjd2000,data,sim);

% 3D plot of the relative position of all ast during coasting
figure('Name','Coasting rel pos ast earth')
plot3(coasting.ast1.rel_pos_ast_earth(:,1),coasting.ast1.rel_pos_ast_earth(:,2),coasting.ast1.rel_pos_ast_earth(:,3),...
    'DisplayName','Ast 1','Color',colors(1,:));
axis equal; hold on;
plot3(coasting.ast2.rel_pos_ast_earth(:,1),coasting.ast2.rel_pos_ast_earth(:,2),coasting.ast2.rel_pos_ast_earth(:,3),...
    'DisplayName','Ast 2','Color',colors(2,:));
plot3(coasting.asta.rel_pos_ast_earth(:,1),coasting.asta.rel_pos_ast_earth(:,2),coasting.asta.rel_pos_ast_earth(:,3),...
    'DisplayName','Ast a','Color',colors(3,:));
plot3(coasting.astb.rel_pos_ast_earth(:,1),coasting.astb.rel_pos_ast_earth(:,2),coasting.astb.rel_pos_ast_earth(:,3),...
    'DisplayName','Ast b','Color',colors(4,:));
plot3(0,0,0,'o','DisplayName','Earth','Color',colors(8,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Coasting - 3D rel pos of Ast wrt EA');

% 3D plot of the relative position of earth wrt all ast, during coasting
figure('Name','Coasting rel pos ast earth')
subplot(2,2,1)
plot3(coasting.ast1.rel_pos_earth_ast(:,1),coasting.ast1.rel_pos_earth_ast(:,2),coasting.ast1.rel_pos_earth_ast(:,3),...
    'DisplayName','Ast 1','Color',colors(1,:));
axis equal; hold on;
plot3(0,0,0,'o','DisplayName','Ast 1','Color',colors(2,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
subplot(2,2,2)
plot3(coasting.ast2.rel_pos_earth_ast(:,1),coasting.ast2.rel_pos_earth_ast(:,2),coasting.ast2.rel_pos_earth_ast(:,3),...
    'DisplayName','Ast 2','Color',colors(1,:));
axis equal; hold on;
plot3(0,0,0,'o','DisplayName','Ast 2','Color',colors(3,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
subplot(2,2,3)
plot3(coasting.asta.rel_pos_earth_ast(:,1),coasting.asta.rel_pos_earth_ast(:,2),coasting.asta.rel_pos_earth_ast(:,3),...
    'DisplayName','Ast a','Color',colors(1,:));
axis equal; hold on;
plot3(0,0,0,'o','DisplayName','Ast 3','Color',colors(4,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
subplot(2,2,4)
plot3(coasting.astb.rel_pos_earth_ast(:,1),coasting.astb.rel_pos_earth_ast(:,2),coasting.astb.rel_pos_earth_ast(:,3),...
    'DisplayName','Ast b','Color',colors(1,:));
axis equal; hold on;
plot3(0,0,0,'o','DisplayName','Ast 4','Color',colors(5,:))
legend('show')
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Coasting - 3D rel pos of EA wrt Ast');

% plot of magnitude of the relative distance ast - earth
figure('Name','magnitude of rel pos')
plot(coasting.ast1.time-sol.departure_mjd2000,coasting.ast1.norm_rel_pos_ast_earth,'DisplayName','Ast 1','Color',colors(1,:))
hold on
plot(coasting.ast2.time-sol.departure_mjd2000,coasting.ast2.norm_rel_pos_ast_earth,'DisplayName','Ast 2','Color',colors(2,:))
plot(coasting.asta.time-sol.departure_mjd2000,coasting.asta.norm_rel_pos_ast_earth,'DisplayName','Ast a','Color',colors(3,:))
plot(coasting.astb.time-sol.departure_mjd2000,coasting.astb.norm_rel_pos_ast_earth,'DisplayName','Ast b','Color',colors(4,:))
legend('show')
xlabel('time from departure [d]'); ylabel('norm(REL DIST ast wrt EA) [AU]'); 

% plot of magnitude of the relative distance earth - ast
figure('Name','magnitude of rel pos')
plot(coasting.ast1.time-sol.departure_mjd2000,coasting.ast1.norm_rel_pos_earth_ast,'DisplayName','Ast 1','Color',colors(1,:))
hold on
plot(coasting.ast2.time-sol.departure_mjd2000,coasting.ast2.norm_rel_pos_earth_ast,'DisplayName','Ast 2','Color',colors(2,:))
plot(coasting.asta.time-sol.departure_mjd2000,coasting.asta.norm_rel_pos_earth_ast,'DisplayName','Ast a','Color',colors(3,:))
plot(coasting.astb.time-sol.departure_mjd2000,coasting.astb.norm_rel_pos_earth_ast,'DisplayName','Ast b','Color',colors(4,:))
legend('show')
xlabel('time from departure [d]'); ylabel('norm(REL DIST EA wrt ast) [AU]'); 

% actual plot of asteroid and earth, heliocentric
figure('Name','ast ea')
plot3(coasting.ast1.r_ast(:,1)./sim.DU,coasting.ast1.r_ast(:,2)./sim.DU,coasting.ast1.r_ast(:,3)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','ast');
axis equal; hold on;
plot3(coasting.ast1.rEA(:,1)./sim.DU,coasting.ast1.rEA(:,2)./sim.DU,coasting.ast1.rEA(:,3)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(2,:),'DisplayName','EA');
plot3(0,0,0,'o','Color',colors(4,:))
legend('show')
view(2)

% cartesian element of the asteroid 1, just to check
figure('Name','ast1 xyz')
plot(coasting.ast1.time, coasting.ast1.r_ast(:,1)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','x');
hold on
plot(coasting.ast1.time, coasting.ast1.r_ast(:,2)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(2,:),'DisplayName','y');
plot(coasting.ast1.time, coasting.ast1.r_ast(:,3)./sim.DU,...
    'DisplayName','Ast 1','Color',colors(3,:),'DisplayName','z');
legend('show')

% orbital parameters plot of 1st ast to check
figure('Name','orbital params')
subplot(3,2,1)
plot(coasting.ast1.time, coasting.ast1.coe_ast(:,1),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','a');
legend('show')
subplot(3,2,2)
plot(coasting.ast1.time, coasting.ast1.coe_ast(:,2),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','e');
legend('show')
subplot(3,2,3)
plot(coasting.ast1.time, coasting.ast1.coe_ast(:,3),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','i');
legend('show')
subplot(3,2,4)
plot(coasting.ast1.time, coasting.ast1.coe_ast(:,4),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','OM');
legend('show')
subplot(3,2,5)
plot(coasting.ast1.time, coasting.ast1.coe_ast(:,5),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','om');
legend('show')
subplot(3,2,6)
plot(coasting.ast1.time, coasting.ast1.coe_ast(:,6),...
    'DisplayName','Ast 1','Color',colors(1,:),'DisplayName','theta');
legend('show')