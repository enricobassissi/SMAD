function [] = plot_mission(sol,colors)
    
    % initialise stuff from solution
    MJD01 = sol.MJD0;
    MJDF1 = MJD01 + sol.TOF1;
    MJD02 = MJDF1 + sol.buffer_time;
    MJDF2 = MJD02 + sol.TOF2;
    
    % position of planet at sol moments
    [kep_EA,ksun] = uplanet(MJD01, 3);
    [rm_EA, ~] = sv_from_coe(kep_EA,ksun);
    [kep_MA1,ksun] = uplanet(MJDF1, 4);
    [rm_MA1, vm_MA1] = sv_from_coe(kep_MA1,ksun);
    [kep_MA2,ksun] = uplanet(MJD02, 4);
    [rm_MA2, ~] = sv_from_coe(kep_MA2,ksun);
    [kep_JU,ksun] = uplanet(MJDF2, 5);
    [rm_JU, ~] = sv_from_coe(kep_JU,ksun);
    
    % Converting mJ2000 in seconds
    t1_sec = MJD01.*60.*60.*24;
    t2_sec = MJDF1.*60.*60.*24;
    t3_sec = MJD02.*60.*60.*24;
    t4_sec = MJDF2.*60.*60.*24;

    ToF12_sec = t2_sec - t1_sec;
    [~,~,~,~,VI12,VF12,~,~] = lambertMR(rm_EA,rm_MA1,ToF12_sec,ksun,0,0,0,0);
    
    ToF34 = t4_sec - t3_sec;
    [~,~,~,~,VI34,~,~,~] = lambertMR(rm_MA2,rm_JU,ToF34,ksun,0,0,0,0);

    % FULL ORBIT IN THE YEAR OF THE ACTUAL ENCOUNTER
    oneyearJU=12*365; % one year of Jupiter in days 
    oneyearMA=687; % one year of Mars
    oneyearEA=365; % one year of Mercury

    % for full orbit propagation
    n=100;
    T_EA=linspace(MJD01,MJD01+oneyearEA,n);
    T_MA1=linspace(MJDF1,MJDF1+oneyearMA,n);
    T_MA2=linspace(MJD02,MJD02+oneyearMA,n);
    T_JU=linspace(MJDF2,MJDF2+oneyearJU,n);


    % INITIAL AND FINAL ORBIT DEFINITION
    % ID JUPITER 5, ID MARS 4, ID Earth 3
    % [kep,ksun] = uplanet(mjd2000, ibody)

    % COMPLETE ORBITS
    R_EA = zeros(n,3); V_EA = zeros(n,3);
    R_MA1 = zeros(n,3); V_MA1 = zeros(n,3);
    R_MA2 = zeros(n,3); V_MA2 = zeros(n,3);
    R_JU = zeros(n,3); V_JU = zeros(n,3);
    for k=1:n
        [kep_EA,ksun] = uplanet(T_EA(k), 3);
        [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);  
        R_EA(k,:)=r_EA;
        V_EA(k,:)=v_EA;
        [kep_MA1,ksun] = uplanet(T_MA1(k), 4);
        [r_MA1, v_MA1] = sv_from_coe(kep_MA1,ksun);  
        R_MA1(k,:)=r_MA1;
        V_MA1(k,:)=v_MA1;
        [kep_MA2,ksun] = uplanet(T_MA2(k), 4);
        [r_MA2, v_MA2] = sv_from_coe(kep_MA2,ksun);  
        R_MA2(k,:)=r_MA2;
        V_MA2(k,:)=v_MA2;
        [kep_JU,ksun] = uplanet(T_JU(k), 5);
        [r_JU, v_JU] = sv_from_coe(kep_JU,ksun);  
        R_JU(k,:)=r_JU;
        V_JU(k,:)=v_JU;
    end

    % PLOT ORBITS AND BEST LAMBERT TRANSFER 
    figure()

    AU = astroConstants(2);

    hold on
    hJU = plot3(R_JU(:,1)./AU,R_JU(:,2)./AU,R_JU(:,3)./AU,'--','Color',colors(5,:));
    hMA1 = plot3(R_MA1(:,1)./AU,R_MA1(:,2)./AU,R_MA1(:,3)./AU,'--','Color',colors(2,:));
    hMA2 = plot3(R_MA2(:,1)./AU,R_MA2(:,2)./AU,R_MA2(:,3)./AU,'-.','Color',colors(2,:));
    hEA = plot3(R_EA(:,1)./AU,R_EA(:,2)./AU,R_EA(:,3)./AU,'--','Color',colors(8,:));

    % first leg
    t012 = t1_sec;
    tf12 = t2_sec;
    y012 = [rm_EA, VI12]'; % velocity from lambert arc transfer orbit injection
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y12] = ode113(@rates, [t012 tf12], y012,options,'sun');
    h1 = plot3( y12(:,1)./AU, y12(:,2)./AU, y12(:,3)./AU,'Color',colors(3,:));

    % coasting
    t023 = t2_sec;
    tf23 = t3_sec;
    y023 = [rm_MA1, vm_MA1]'; % after the application of 2nd dv of first leg, rendezvous with mars, same velocity as mars at that moment
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y23] = ode113(@rates, [t023 tf23], y023,options,'sun');
    h2 = plot3( y23(:,1)./AU, y23(:,2)./AU, y23(:,3)./AU,'*','Color',colors(3,:),'Markersize',3);

    % second leg
    t034 = t3_sec;
    tf34 = t4_sec;
    y034 = [rm_MA2, VI34]';
    options = odeset ('RelTol', 1e-13, 'AbsTol', 1e-14); 
    [~,y34] = ode113(@rates, [t034 tf34], y034,options,'sun');
    h3 = plot3( y34(:,1)./AU, y34(:,2)./AU, y34(:,3)./AU,'-.','Color',colors(3,:));

    h4 = plot3(0,0,0,'*','Color',colors(4,:));
    plot3(rm_EA(1)./AU,rm_EA(2)./AU,rm_EA(3)./AU,'o','Color',colors(8,:),'MarkerSize',4);
    plot3(rm_MA1(1)./AU,rm_MA1(2)./AU,rm_MA1(3)./AU,'o','Color',colors(2,:),'MarkerSize',4);
    plot3(rm_MA2(1)./AU,rm_MA2(2)./AU,rm_MA2(3)./AU,'o','Color',colors(2,:),'MarkerSize',4);
    plot3(rm_JU(1)./AU,rm_JU(2)./AU,rm_JU(3)./AU,'o','Color',colors(5,:));

    axis equal; grid on
    legend ([hEA,hMA1,hJU,h1,h2,h3,h4],'Full Earth','Full Mars','Full Jupiter',...
            'Mission Leg 1','Coasting','Mission Leg 2','Sun','Location','southeast')
    xlabel('AU')
    ylabel('AU')
    zlabel('AU')

end