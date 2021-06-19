function [actual_elements] = uNEO2(mjd2000,name,data)

%{
 This function interpolate the ephemerides contained in data.y_interp_ft.

 data.y_interp_ft has been calculated outside of this code to speed up the
 process and it's derived from a fourier interpolation of the real DE431 JPL
 Horizons Ephemerides from before MORPHEUS mission starting date (2020) to
 after the mission end date (2050), such that each possible date of
 analysis fall inside the boundaries.

 For each asteroid, a cell of [200x6] elements is built, with its related time
 vector of evaluation data.t_vector [1x200].

 data contains:
 asteroid_names:            the list of our selected asteroid [9x1 string]
 data.PermutationMatrix:    the unique permutation matrix of 4 asteroids among 9 
                            elements [3024x4 string] 
 data.y_interp_ft:          cells containing the [200x6] interpolating parameters for 
                            all the considered asteroids [9x1 cell]
 data.t_vector:             time vector of evaluation of the above points [1x200 double]

 The data.y_interp_ft from struct data contains:
 [a,e,i,OM,om,theta] in [km,-,rad,rad,rad,wrapped rad]
%}

    % Find which asteroid we want
    idx_asteroid_in_the_list = name == data.asteroid_names;

    % Extrapolate its orbital elements matrix
    FT_COEFF = data.y_interp_ft{idx_asteroid_in_the_list};

    % Interpolate the actual time we want with a spline between two consecutive lookup-table points
    actual_elements = interp1(data.t_vector, FT_COEFF, mjd2000, 'spline'); %interpolate data with cubic line

end