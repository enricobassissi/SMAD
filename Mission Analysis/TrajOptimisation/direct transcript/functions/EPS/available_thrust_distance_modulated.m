function [T_max_modulated] = available_thrust_distance_modulated(T_max, RR)
%	available_thrust_distance_modulated - Modulate the maximum thrust with
%	the distance
%
% PROTOTYPE:
%  [T_max_modulated] = available_thrust_distance_modulated(T_max, RR)
%
% DESCRIPTION:
%   Compute the thust available for the mission
%
%  INPUT :
%   T_max:      [1] maximum thrust
%   RR:     [Nx3] Array of vectors position in heliocentric cartesian frame [AU]
%
%  OUTPUT:
%   T_max_modulated:                       [Nx1] MAx thrust modulated [mN]
%
%  FUNCTIONS CALLED:
%   solar_irradiance_cooler
%
% AUTHOR:
%   Marco Elia
%
% PREVIOUS VERSION:
%   \\
%
% CHANGELOG:
%   23/06/2021, Marco Elia
%
% -------------------------------------------------------------------------
% REFERENCES:
%
d_sun = vecnorm(RR,2,2);
eta_space = solar_irradiance_cooler(d_sun, 1, 1);
eta_space = eta_space/max(eta_space);
T_max_modulated = T_max*eta_space;
end






