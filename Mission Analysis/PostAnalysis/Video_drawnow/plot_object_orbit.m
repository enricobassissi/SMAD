function [R, h] = plot_object_orbit(mjd2000,obj_id,period,sim,data,n,colors,color_id,varargin)
    
    AU = sim.DU;
    muSun = sim.mu_dim;
    
    switch lower(obj_id)
        case 'earth'
            % Earth Full Orbit
            disp_name = 'Earth';
%             oneyear = 365; %d
            T = linspace(mjd2000-period/2,mjd2000+period/2,n);
            for k=1:n
                [kep,ksun] = uplanet(T(k), 3);
                [r, ~] = sv_from_coe(kep,ksun);  
                R(k,:)=r/AU; % it's in km it becomes AU, to be plotted
            end
        case 'mars'
            % Mars Full Orbit
            disp_name = 'Mars';
%             oneyear = 687; %d
            T = linspace(mjd2000-period/2,mjd2000+period/2,n);
            for k=1:n
                [kep,ksun] = uplanet(T(k), 4);
                [r, ~] = sv_from_coe(kep,ksun);  
                R(k,:)=r/AU; % it's in km it becomes AU, to be plotted
            end
        case 'asteroid'
            % Ast Full Orbit
            disp_name = varargin{1};
%             py_ast = py.list(disp_name);
%             % Query JPL SBDB for the bodies in the risk list
%             py_dict_ast=py.neo_api_function.get_dict(py_ast);
%             % Selected Asteroids Characteristics Cell
%             [selected_asteroids_orbital_elements_and_sigma, ~] = ...
%                 get_orbital_elements_and_sigma(ast_name,py_dict_ast);
%             % Period of the Ast
%             a_ast = selected_asteroids_orbital_elements_and_sigma{1}(1,1)*AU; % km
%             T_ast = 2*pi*sqrt(a_ast^3/muSun); % s
%             T_ast_days = T_ast/(3600*24); % days
            T = linspace(mjd2000-period/2,mjd2000+period/2,n);
            disp_name = string(disp_name);
            for k=1:n
                [kep_ast] = uNEO3(T(k),disp_name,data); % [km,-,rad,rad,rad,wrapped rad]
                [r, ~] = sv_from_coe(kep_ast,muSun); % km, km/s
                R(k,:) = r/AU;
            end
    end

    h = plot3(R(:,1),R(:,2),R(:,3),'--','Color',colors(color_id,:),...
        'LineWidth',1); %,'DisplayName', disp_name
    
end