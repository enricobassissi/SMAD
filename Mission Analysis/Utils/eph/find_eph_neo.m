function [y_interp_ft, time_vector] = find_eph_neo(asteroid_list)
    
    % output:
    % [a,e,i,OM,om,theta] in [km,-,rad,rad,rad,wrapped rad]
    % time_vector [days]
    
    AU = 149597870.691; % km From DE405
    
    PointOfView = 'Sun';
    epoch_start = '2020-01-01';
    epoch_stop = '2050-01-01';
    step = '1d';
    type_elements = 'elements';
    
    y_interp_ft = cell(length(asteroid_list),1);
    for name = 1:length(asteroid_list)
        % data extraction section
        py_data = py.neo_api_function.get_horizons_ephemerides(py.str(asteroid_list(name)),...
            py.str(PointOfView),py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
        horizons_data = double(py_data); % [a] in AU; [i,OM,om,theta] in deg
        
        % unwrapping of theta for better interpolation
        horizons_data(:,6) = unwrap(deg2rad(horizons_data(:,6))); % rad
        
        % up to this pint we have [a,e,i,OM,om,theta] in [AU,-,deg,deg,deg,unwrapped rad]
        
        % Fourier interpolation
        N = 200; % number of query point, interpolation points
        y_interp_ft{name,1} = interpft(horizons_data,N,1); % 1 by column, 2 by rows
        
        % wrap again to 2pi the true anomaly
        y_interp_ft{name,1}(:,6) = wrapping_to_2pi(y_interp_ft{name,1}(:,6)); % wrapped rad
        % change output dimensions
        y_interp_ft{name,1}(:,1) = y_interp_ft{name,1}(:,1)*AU; % AU -> km
        y_interp_ft{name,1}(:,3:5) = deg2rad(y_interp_ft{name,1}(:,3:5)); % deg -> rad
        % true anomaly y_interp(6) is already in rad
    end
    
    time_vector = linspace(pystr2mjd2000(epoch_start),pystr2mjd2000(epoch_stop),N); %days

end