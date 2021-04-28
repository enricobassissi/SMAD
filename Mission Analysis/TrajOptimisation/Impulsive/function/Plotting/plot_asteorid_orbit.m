function [R_AST, h] = plot_asteorid_orbit(mjd2000,years,ast_name,colors,color_id)

    epoch_start = mjd20002pystr(mjd2000); 
    epoch_stop = mjd20002pystr(mjd2000+years*365);
    step = '10d'; 
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(ast_name),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data = double(py_data);
    R_AST = horizons_data(:, 1:3); %AU
    
    h = plot3(R_AST(:,1),R_AST(:,2),R_AST(:,3),'--','Color',colors(color_id,:),...
        'DisplayName',strcat('Full ',ast_name));

end