function [R_AST] = coasting_asteroids_2(mjd2000,mjd2000_final,ast_name)

%     epoch_start = mjd20002pystr(mjd2000-round(T_ast_days*fraction_of_the_orbit/2)); 
%     epoch_stop = mjd20002pystr(mjd2000+round(T_ast_days*fraction_of_the_orbit/2));
    epoch_start = mjd20002pystr(mjd2000); 
    epoch_stop = mjd20002pystr(mjd2000_final);
    step = '1d'; 
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(ast_name),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data = double(py_data);
    R_AST = horizons_data(:, 1:3); %AU
%     
%     h = plot3(R_AST(:,1),R_AST(:,2),R_AST(:,3),'--','Color',colors(color_id,:),'LineWidth',0.8);
% %         'DisplayName',strcat('Full ',ast_name)
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';

end