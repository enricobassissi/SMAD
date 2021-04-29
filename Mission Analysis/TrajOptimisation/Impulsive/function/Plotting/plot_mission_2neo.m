function [] = plot_mission_2neo(sol,colors,ast_name)
    
    AU = astroConstants(2);

    % initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDF1 = MJD01 + sol.TOF1;
    MJD02 = MJDF1 + sol.buffer_time;
    MJDF2 = MJD02 + sol.TOF2;
    
    ast1 = ast_name(1);
    ast2 = ast_name(2);
    
    % position of planet at sol moments
    [kep_EA, muSun] = uplanet(MJD01, 3);
    [rm_EA, ~] = sv_from_coe(kep_EA,muSun); %km, km/s
    [r_ast1_1,v_ast1_1] = uNEO(MJDF1, ast1); %km, km/s
    [r_ast1_2,~] = uNEO(MJD02, ast1);
    [r_ast2,~] = uNEO(MJDF2, ast2);
    
    % Converting mJ2000 in seconds
    t1_sec = MJD01.*60.*60.*24;
    t2_sec = MJDF1.*60.*60.*24;
    t3_sec = MJD02.*60.*60.*24;
    t4_sec = MJDF2.*60.*60.*24;
    
    % Lamberts
    ToF12_sec = t2_sec - t1_sec;
    [~,~,~,~,VI12,~,~,~] = lambertMR(rm_EA,r_ast1_1,ToF12_sec,muSun,0,0,0,0);
    
    ToF34 = t4_sec - t3_sec;
    [~,~,~,~,VI34,~,~,~] = lambertMR(r_ast1_2,r_ast2,ToF34,muSun,0,0,0,0);

    % FULL ORBIT IN THE YEAR OF THE ACTUAL ENCOUNTER
    % Asteroids
    epoch_start = mjd20002pystr(MJDF1); 
    epoch_stop = mjd20002pystr(MJDF1+6*365);
    step = '10d'; % to produce 2 eph element: eph for day (MJDF1-1) at midnight and one for day (MJDF1) at midnight
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data1 = py.neo_api_function.get_horizons_ephemerides(py.str(ast1),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data1 = double(py_data1);
    R_AST1 = horizons_data1(:, 1:3); %AU
    
    epoch_start = mjd20002pystr(MJDF2); 
    epoch_stop = mjd20002pystr(MJDF2+6*365);
    step = '10d'; % to produce 2 eph element: eph for day (MJDF1-1) at midnight and one for day (MJDF1) at midnight
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data2 = py.neo_api_function.get_horizons_ephemerides(py.str(ast2),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data2 = double(py_data2);
    R_AST2 = horizons_data2(:, 1:3); %AU
    
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

    % PLOT ORBITS AND BEST LAMBERT TRANSFER 
    figure('Name','Mission Orbits and Phases')

    hold on
    plot3(R_AST2(:,1),R_AST2(:,2),R_AST2(:,3),'--','Color',colors(5,:),...
        'DisplayName',strcat('Full ',ast2));
    plot3(R_AST1(:,1),R_AST1(:,2),R_AST1(:,3),'--','Color',colors(2,:),...
        'DisplayName',strcat('Full ',ast1));
    plot3(R_EA(:,1),R_EA(:,2),R_EA(:,3),'--','Color',colors(8,:),...
        'DisplayName','Full Earth');

    % first leg
    t012 = t1_sec;
    tf12 = t2_sec;
    y012 = [rm_EA, VI12]'; %km, km/s; velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y12] = ode113(@rates, [t012 tf12], y012,options,'sun');
    plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(3,:),...
        'DisplayName','First Leg');

    % coasting
    t023 = t2_sec;
    tf23 = t3_sec;
    y023 = [r_ast1_1, v_ast1_1]'; %km, km/s; after the application of 2nd dv of first leg, rendezvous with mars, same velocity as mars at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y23] = ode113(@rates, [t023 tf23], y023,options,'sun');
    plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'*','Color',colors(3,:),'Markersize',3,...
        'DisplayName','Coasting');

    % second leg
    t034 = t3_sec;
    tf34 = t4_sec;
    y034 = [r_ast1_2, VI34]'; %km, km/s
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y34] = ode113(@rates, [t034 tf34], y034,options,'sun');
    plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'-.','Color',colors(3,:),...
        'DisplayName','Leg 2');

    plot3(0,0,0,'*','Color',colors(4,:),'DisplayName','Sun');
    
    legend('show','Location','southeastoutside')
    
    hp1 = plot3(rm_EA(1)./AU,rm_EA(2)./AU,rm_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',4);
    hp1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp2 = plot3(r_ast1_1(1)./AU,r_ast1_1(2)./AU,r_ast1_1(3)./AU,'o','Color',colors(2,:),'MarkerSize',4);
    hp2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp3 = plot3(r_ast1_2(1)./AU,r_ast1_2(2)./AU,r_ast1_2(3)./AU,'o','Color',colors(2,:),'MarkerSize',4);
    hp3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp4 = plot3(r_ast2(1)./AU,r_ast2(2)./AU,r_ast2(3)./AU,'o','Color',colors(5,:));
    hp4.Annotation.LegendInformation.IconDisplayStyle = 'off';

    axis equal; grid on
    xlabel('AU')
    ylabel('AU')
    zlabel('AU')

end