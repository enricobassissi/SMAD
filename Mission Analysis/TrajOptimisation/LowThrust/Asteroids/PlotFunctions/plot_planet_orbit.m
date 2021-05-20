function [R, h] = plot_planet_orbit(mjd2000,PL_ID,colors,color_id)
    
    AU = astroConstants(2);
    n=100;
    
    switch PL_ID 
        case 3
            % Earth Full Orbit
            disp_name = 'Earth';
            oneyear = 365; %d
            T = linspace(mjd2000,mjd2000+oneyear,n);
            R = zeros(n,3); V = zeros(n,3);
        case 4
            % Mars Full Orbit
            disp_name = 'Mars';
            oneyear = 687; %d
            T = linspace(mjd2000,mjd2000+oneyear,n);
            R = zeros(n,3); V = zeros(n,3);
    end

    for k=1:n
        [kep,ksun] = uplanet(T(k), PL_ID);
        [r, v] = sv_from_coe(kep,ksun);  
        R(k,:)=r/AU; % it's in km it becomes AU, to be plotted
        V(k,:)=v;
    end

    h = plot3(R(:,1),R(:,2),R(:,3),'--','Color',colors(color_id,:),...
        'LineWidth',1,'DisplayName', disp_name);
    
end