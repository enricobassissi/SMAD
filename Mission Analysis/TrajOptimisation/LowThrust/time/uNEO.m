function [r,v] = uNEO(mjd2000,name)

    AU = astroConstants(2); %km
    day = 86400; %s

    epoch_start = mjd20002pystr(mjd2000-1); %just to produce the ephemerides
    epoch_stop = mjd20002pystr(mjd2000);
    step = '1d'; % to produce 2 eph element: eph for day (MJDF1-1) at midnight and one for day (MJDF1) at midnight
    type_elements = 'Vectors';
    PointOfView = 'Sun';
    asteroid_name = name; % string name
    py_data = py.neo_api_function.get_horizons_ephemerides(py.str(asteroid_name),py.str(PointOfView),...
                      py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
    horizons_data = double(py_data);
    
    % We take the last of the generated ephemerides that relates to the actual MJDF1
    r = horizons_data(end, 1:3)'*AU; %it's in AU -> and it becomes km
    v = horizons_data(end, 4:6)'*AU/day; %it's in AU/day -> and it becomes km/s
    
end