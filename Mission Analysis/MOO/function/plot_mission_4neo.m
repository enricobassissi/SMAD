function [] = plot_mission_4neo(sol,colors,asteroid_names_sequence)
    
    AU = astroConstants(2);

    % initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDF1 = MJD01 + sol.TOF1;
    MJD02 = MJDF1 + sol.buffer_time1;
    MJDF2 = MJD02 + sol.TOF2;
    MJD03 = MJDF2 + sol.buffer_time2;
    MJDF3 = MJD03 + sol.TOF3;
    MJD04 = MJDF3 + sol.buffer_time3;
    MJDF4 = MJD04 + sol.TOF4;
    
    ast1 = asteroid_names_sequence(1);
    ast2 = asteroid_names_sequence(2);
    ast3 = asteroid_names_sequence(3);
    ast4 = asteroid_names_sequence(4);
    
    % Position of planet at solutions moments
    [kep_EA, muSun] = uplanet(MJD01, 3); % Earth Departure
    [rm_EA, ~] = sv_from_coe(kep_EA,muSun); %km, km/s
    [r_ast1_1,v_ast1_1] = uNEO(MJDF1, ast1); %km, km/s, Ast 1 Arrival
    [r_ast1_2,~] = uNEO(MJD02, ast1); % Ast 1 Departure after Coasting
    [r_ast2_1,v_ast2_1] = uNEO(MJDF2, ast2); % Ast 2 Arrival
    [r_ast2_2,~] = uNEO(MJD03, ast2); % Ast 2 Departure after Coasting
    [r_ast3_1,v_ast3_1] = uNEO(MJDF3, ast3); % Ast 3 Arrival
    [r_ast3_2,~] = uNEO(MJD04, ast3); % Ast 3 Departure after Coasting
    [r_ast4_1,~] = uNEO(MJDF4, ast4); % Ast 4 Arrival
    
    % Converting mJD2000 in seconds
    t1_sec = MJD01*60*60*24;
    t2_sec = MJDF1*60*60*24;
    t3_sec = MJD02*60*60*24;
    t4_sec = MJDF2*60*60*24;
    t5_sec = MJD03*60*60*24;
    t6_sec = MJDF3*60*60*24;
    t7_sec = MJD04*60*60*24;
    t8_sec = MJDF4*60*60*24;
    
    % Lamberts
    ToF12_sec = t2_sec - t1_sec;
    [~,~,~,~,VI12,~,~,~] = lambertMR(rm_EA,r_ast1_1,ToF12_sec,muSun,0,0,0,0);
    ToF34_sec = t4_sec - t3_sec;
    [~,~,~,~,VI34,~,~,~] = lambertMR(r_ast1_2,r_ast2_1,ToF34_sec,muSun,0,0,0,0);
    ToF56_sec = t6_sec - t5_sec;
    [~,~,~,~,VI56,~,~,~] = lambertMR(r_ast2_2,r_ast3_1,ToF56_sec,muSun,0,0,0,0);
    ToF78_sec = t8_sec - t7_sec;
    [~,~,~,~,VI78,~,~,~] = lambertMR(r_ast3_2,r_ast4_1,ToF78_sec,muSun,0,0,0,0);

    % FULL ORBIT IN THE YEAR OF THE ACTUAL ENCOUNTER
    % Earth Full Orbit
    n=100;
    oneyearEA = 365; %d
    T_EA = linspace(MJD01,MJD01+oneyearEA,n);
    R_EA = zeros(n,3); V_EA = zeros(n,3);
    for k=1:n
        [kep_EA,muSun] = uplanet(T_EA(k), 3);
        [r_EA, v_EA] = sv_from_coe(kep_EA,muSun);  
        R_EA(k,:)=r_EA/AU; % it's in km it becomes already AU, to be plotted
        V_EA(k,:)=v_EA;
    end
    
    % Asteroids
    epoch_start = mjd20002pystr(MJDF1); 
    epoch_stop = mjd20002pystr(MJDF1+6*365);
    step = '10d'; 
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data1 = py.neo_api_function.get_horizons_ephemerides(py.str(ast1),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data1 = double(py_data1);
    R_AST1 = horizons_data1(:, 1:3); %AU
    
    epoch_start = mjd20002pystr(MJDF2); 
    epoch_stop = mjd20002pystr(MJDF2+6*365);
    step = '10d';
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data2 = py.neo_api_function.get_horizons_ephemerides(py.str(ast2),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data2 = double(py_data2);
    R_AST2 = horizons_data2(:, 1:3); %AU
    
    epoch_start = mjd20002pystr(MJDF3); 
    epoch_stop = mjd20002pystr(MJDF3+6*365);
    step = '10d'; 
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data3 = py.neo_api_function.get_horizons_ephemerides(py.str(ast3),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data3 = double(py_data3);
    R_AST3 = horizons_data3(:, 1:3); %AU
    
    epoch_start = mjd20002pystr(MJDF4); 
    epoch_stop = mjd20002pystr(MJDF4+6*365);
    step = '10d'; 
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data4 = py.neo_api_function.get_horizons_ephemerides(py.str(ast4),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data4 = double(py_data4);
    R_AST4 = horizons_data4(:, 1:3); %AU

    % PLOT ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')
    hold on
    plot3(R_EA(:,1),R_EA(:,2),R_EA(:,3),'--','Color',colors(8,:),...
        'DisplayName','Full Earth');
    plot3(R_AST1(:,1),R_AST1(:,2),R_AST1(:,3),'--','Color',colors(2,:),...
        'DisplayName',strcat('Full ',ast1));
    plot3(R_AST2(:,1),R_AST2(:,2),R_AST2(:,3),'--','Color',colors(3,:),...
        'DisplayName',strcat('Full ',ast2));
    plot3(R_AST3(:,1),R_AST3(:,2),R_AST3(:,3),'--','Color',colors(4,:),...
        'DisplayName',strcat('Full ',ast3));
    plot3(R_AST4(:,1),R_AST4(:,2),R_AST4(:,3),'--','Color',colors(5,:),...
        'DisplayName',strcat('Full ',ast4));
    
    % Mission Arcs
    % First leg: Earth -> Ast 1
    t012 = t1_sec;
    tf12 = t2_sec;
    y012 = [rm_EA, VI12]'; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y12] = ode113(@rates, [t012 tf12], y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','First Leg');

    % Coasting 1 on Ast 1
    t023 = t2_sec;
    tf23 = t3_sec;
    y023 = [r_ast1_1, v_ast1_1]'; %km, km/s; after the application of 2nd dv of its leg, rendezvous with the asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y23] = ode113(@rates, [t023 tf23], y023,options,'sun');
    plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'*','Color',colors(1,:),'Markersize',3,...
        'DisplayName','Coasting 1');

    % Second leg: Ast 1 -> Ast 2
    t034 = t3_sec;
    tf34 = t4_sec;
    y034 = [r_ast1_2, VI34]'; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y34] = ode113(@rates, [t034 tf34], y034,options,'sun');
    plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Leg 2');
    
    % Coasting 2 on Ast 2
    t045 = t4_sec;
    tf45 = t5_sec;
    y045 = [r_ast2_1, v_ast2_1]'; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y45] = ode113(@rates, [t045 tf45], y045,options,'sun');
    plot3( y45(:,1)./AU, y45(:,2)./AU, y45(:,3)./AU,'*','Color',colors(1,:),'Markersize',3,...
        'DisplayName','Coasting 2');

    % Third leg: Ast 2 -> Ast 3
    t056 = t5_sec;
    tf56 = t6_sec;
    y056 = [r_ast2_2, VI56]'; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y56] = ode113(@rates, [t056 tf56], y056,options,'sun');
    plot3( y56(:,1)./AU, y56(:,2)./AU, y56(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Leg 3');
    
    % Coasting 3 on Ast 3
    t067 = t6_sec;
    tf67 = t7_sec;
    y067 = [r_ast3_1, v_ast3_1]'; %km, km/s; after the application of 2nd dv of its transfer leg, rendezvous with asteroid, same velocity as the asteroid at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y67] = ode113(@rates, [t067 tf67], y067,options,'sun');
    plot3( y67(:,1)./AU, y67(:,2)./AU, y67(:,3)./AU,'*','Color',colors(1,:),'Markersize',3,...
        'DisplayName','Coasting 3');

    % Fourth leg: Ast 3 -> Ast 4
    t078 = t7_sec;
    tf78 = t8_sec;
    y078 = [r_ast3_2, VI78]'; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y78] = ode113(@rates, [t078 tf78], y078,options,'sun');
    plot3( y78(:,1)./AU, y78(:,2)./AU, y78(:,3)./AU,'Color',colors(1,:),...
        'DisplayName','Leg 4');
    
    % Sun Yellow Asterisk
    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(rm_EA(1)./AU,rm_EA(2)./AU,rm_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',4);
    hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp2 = plot3(r_ast1_1(1)./AU,r_ast1_1(2)./AU,r_ast1_1(3)./AU,'^','Color',colors(2,:),'MarkerSize',4);
    hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp3 = plot3(r_ast1_2(1)./AU,r_ast1_2(2)./AU,r_ast1_2(3)./AU,'o','Color',colors(2,:),'MarkerSize',4);
    hp3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp4 = plot3(r_ast2_1(1)./AU,r_ast2_1(2)./AU,r_ast2_1(3)./AU,'^','Color',colors(3,:),'MarkerSize',4);
    hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp5 = plot3(r_ast2_2(1)./AU,r_ast2_2(2)./AU,r_ast2_2(3)./AU,'o','Color',colors(3,:),'MarkerSize',4);
    hp5.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp6 = plot3(r_ast3_1(1)./AU,r_ast3_1(2)./AU,r_ast3_1(3)./AU,'^','Color',colors(4,:),'MarkerSize',4);
    hp6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp7 = plot3(r_ast3_2(1)./AU,r_ast3_2(2)./AU,r_ast3_2(3)./AU,'o','Color',colors(4,:),'MarkerSize',4);
    hp7.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp8 = plot3(r_ast4_1(1)./AU,r_ast4_1(2)./AU,r_ast4_1(3)./AU,'^','Color',colors(5,:),'MarkerSize',4);
    hp8.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
    xlabel('AU')
    ylabel('AU')
    zlabel('AU')

end