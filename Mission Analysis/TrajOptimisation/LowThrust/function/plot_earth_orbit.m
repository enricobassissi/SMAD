function [R_EA, h] = plot_earth_orbit(mjd2000,colors,color_id)
    
    AU = astroConstants(2);
    
    % Earth Full Orbit
    n=100;
    oneyearEA = 365; %d
    T_EA = linspace(mjd2000,mjd2000+oneyearEA,n);
    R_EA = zeros(n,3); V_EA = zeros(n,3);
    for k=1:n
        [kep_EA,ksun] = uplanet(T_EA(k), 3);
        [r_EA, v_EA] = sv_from_coe(kep_EA,ksun);  
        R_EA(k,:)=r_EA/AU; % it's in km it becomes AU, to be plotted
        V_EA(k,:)=v_EA;
    end

    h = plot3(R_EA(:,1),R_EA(:,2),R_EA(:,3),'--','Color',colors(color_id,:),...
        'DisplayName','Full Earth');
    
end