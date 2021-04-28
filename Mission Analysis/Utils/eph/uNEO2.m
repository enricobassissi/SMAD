function [actual_elements] = uNEO2(mjd2000,name,data)

    % input from data:
    % [a,e,i,OM,om,theta] in [km,-,rad,rad,rad,wrapped rad]
    
    idx_asteroid_in_the_list = name == data.asteroid_names;

    FT_COEFF = data.y_interp_ft{idx_asteroid_in_the_list};

    actual_elements = interp1(data.t_vector, FT_COEFF, mjd2000, 'spline'); %interpolate data with cubic line

end