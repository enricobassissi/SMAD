function [actual_elements] = uNEO3(mjd2000,name,data)

%{
 This function interpolate the ephemerides contained in data.y_interp_ft.

 data.y_interp_ft has been calculated outside of this code to speed up the
 process and it's derived from a fourier interpolation of the real DE431 JPL
 Horizons Ephemerides from before MORPHEUS mission starting date (2020) to
 after the mission end date (2050), such that each possible date of
 analysis fall inside the boundaries.

 N: size of the interpolated data. May differ between a e i OM om theta
    in this case is N=400 for the first 5 elements and N=1000 for the theta
 n: size of the asteroid pool selection. 
 
 For each asteroid, a cell of [Nx6] elements is built, with its related time
 vector of evaluation data.t_vector [1xN].

 data structure contains:
 data.asteroid_names:       the list of our selected asteroid [nx1 string]
 data.PermutationMatrix:    the unique permutation matrix of p asteroids among n 
                            elements [(n!/(n-p)!)x4 string] or depend on
                            the local pruning adopted
 data.y_interp_ft_1_5:      cells containing the [400x6] interpolating parameters for 
                            all the considered asteroids [nx1 cell]
 data.y_interp_ft_6:        cells containing [1000x1] theta element
                            interpolated through fourier
 data.t_vector_1_5:         time vector of evaluation of the above points [1x400 double]
 data.t_vector_6:           time vector of reference for the theta
                            interpolation [1x1000]

 The data.y_interp_ft from struct data contains:
 [a,e,i,OM,om,theta] in [km,-,rad,rad,rad,unwrapped rad]
%}

    % Find which asteroid we want
    idx_asteroid_in_the_list = name == data.asteroid_names;

    % all this because now it's 1-5 discretised with 400 points and the 6th
    % -> theta is with N = 1000
    % Extrapolate its orbital elements matrix
    % [a e i OM om theta]
    FT_COEFF_1_5 = data.y_interp_ft_1_5{idx_asteroid_in_the_list};
    FT_COEFF_6 = data.y_interp_ft_6{idx_asteroid_in_the_list};

    % Interpolate the actual time we want with a spline between two consecutive lookup-table points
    % a is in [km], e stays how the fuck it is, i,OM,om are in [rad],
    % theta is in [rad] too and needs to be wrapped to 2pi
    actual_elements(1,1:5) = interp1(data.t_vector_1_5, FT_COEFF_1_5(:,1:5), mjd2000, 'spline'); %interpolate data with cubic line
    actual_elements(1,6) = wrapTo2Pi(interp1(data.t_vector_6, FT_COEFF_6, mjd2000, 'spline')); %interpolate data with cubic line
end