function [actual_elements] = uNEO2(mjd2000,name,data)

    idx_asteroid_in_the_list = name == data.asteroid_names;

    FT_COEFF = data.y_interp_ft{idx_asteroid_in_the_list};

    actual_elements = interp1(data.t_vector, FT_COEFF, mjd2000, 'spline'); %interpolate data with cubic line
    actual_elements(6) = wrapping_to_2pi(actual_elements(6));

end