function [] = plot_orbits_asteroids(selected_asteroids_names,colors)
    
    figure('Name','All Chosen Asteroids Orbits')
    % Earth 1 year ephemerides
    epoch_start = '2021-01-01';
    epoch_stop = '2022-01-01';
    step = '5d';
    type_elements = 'Vectors';

    py_data_Earth = py.neo_api_function.get_earth_ephemerides(py.str(epoch_start),...
                    py.str(epoch_stop),py.str(step),py.str(type_elements));
    data_Earth = double(py_data_Earth);
    h_EA = plot3(data_Earth(:,1),data_Earth(:,2),data_Earth(:,3),...
               'LineWidth',2.2,"Color",colors(1,:),'DisplayName','Earth');
    % Asteroids
    PointOfView = 'Sun';
    epoch_start = '2021-01-01';
    epoch_stop = '2026-01-01';
    step = '5d';
    type_elements = 'Vectors';

    horizons_data = cell(length(selected_asteroids_names),1);
    for name = 1:length(selected_asteroids_names)
        % data extraction section
        py_data = py.neo_api_function.get_horizons_ephemerides(py.str(selected_asteroids_names(name)),py.str(PointOfView),...
                          py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
        horizons_data{name} = double(py_data);

        % plotting section
%         h_fig = figure();
        
        hold on
        h_asteroid = plot3(horizons_data{name}(:,1),horizons_data{name}(:,2),horizons_data{name}(:,3),...
                     'LineWidth',2,"Color",colors(name+1,:),'DisplayName',selected_asteroids_names(name));
        axis equal; grid on;
        xlabel('x [AU]'); ylabel('y [AU]');  zlabel('z [AU]'); 
%         legend([h_EA, h_asteroid],'Earth',selected_asteroids_names(name),'Location',"northeast");
%         saveas(h_fig,sprintf('./Figures/Orbits/orb%d.png',name));
        %print(h_fig,sprintf('./Figures/Orbits/orb%d.pdf',name),'-dpdf','-bestfit'); 
        %exportgraphics(gca,sprintf('./Figures/Orbits/orb%d.png',name),'ContentType','image');
    end
%     py_data = py.neo_api_function.get_horizons_ephemerides(py.str(selected_asteroids_names(name)),py.str(PointOfView),...
%                          py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
%     horizons_data = double(py_data);
    hold off
    legend show
    
end