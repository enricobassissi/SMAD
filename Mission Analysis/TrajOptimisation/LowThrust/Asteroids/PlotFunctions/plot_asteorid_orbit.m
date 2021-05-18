function [R_AST, h] = plot_asteorid_orbit(mjd2000,fraction_of_the_orbit,ast_name,colors,color_id)
    
    AU = astroConstants(2); %km
    muSun = astroConstants(4); % km^3/s^2
    
    %% get semi major axis of asteroid
    py_ast = py.list({ast_name});
    % Query JPL SBDB for the bodies in the risk list
    py_dict_ast=py.neo_api_function.get_dict(py_ast);
    % Selected Asteroids Characteristics Cell
    [selected_asteroids_orbital_elements_and_sigma, ~] = ...
        get_orbital_elements_and_sigma(ast_name,py_dict_ast);
    % Period of the Ast
    a_ast = selected_asteroids_orbital_elements_and_sigma{1}(1,1)*AU; % km
    T_ast = 2*pi*sqrt(a_ast^3/muSun); % s
    T_ast_days = T_ast/(3600*24);
    
    epoch_start = mjd20002pystr(mjd2000-round(T_ast_days*fraction_of_the_orbit/2)); 
    epoch_stop = mjd20002pystr(mjd2000+round(T_ast_days*fraction_of_the_orbit/2));
    step = '5d'; 
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(ast_name),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data = double(py_data);
    R_AST = horizons_data(:, 1:3); %AU
    
    h = plot3(R_AST(:,1),R_AST(:,2),R_AST(:,3),'--','Color',colors(color_id,:),'LineWidth',0.8);
%         'DisplayName',strcat('Full ',ast_name)
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

end