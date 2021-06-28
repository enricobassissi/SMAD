function [T_max_modulated] = available_thrust_distance_time_angle_modulated(T_max, t, RR, TT)
%	available_thrust_distance_time_angle_modulated - Modulate the maximum thrust with
%	the distance, degradation and angle
%
% PROTOTYPE:
%  [T_max_modulated] = available_thrust_distance_time_angle_modulated(t, RR, TT, data)
%
% DESCRIPTION:
%   Compute the thrust available for the
%	propulsion system during the mission
%
%  INPUT :
%   T_max:      [1] maximum thrust [xN]
%   t:      [Nx1] Array of times [years]
%   RR:     [Nx3] Array of vectors position in heliocentric cartesian frame [AU]
%   TT:     [Nx3] Array of vectors thrust (needed just for the orientation)]
%
%  OUTPUT:
%   T_max_modulated:                       [Nx1] MAx thrust modulated [xN]
%
%  FUNCTIONS CALLED:
%   solar_irradiance_cooler
% AUTHOR:
%   Marco Elia
%
% PREVIOUS VERSION:
%   \\
%
% CHANGELOG:
%   22/06/2021, Marco Elia
%
% -------------------------------------------------------------------------
% REFERENCES:
D = 0.01674;
%distance from the sun
d_sun = vecnorm(RR,2,2);
eta_space = solar_irradiance_cooler(d_sun, 1, 1);

% angle of panels respect to sun
a_sun = acos(dot(RR,TT,2)./(vecnorm(TT,2,2).*d_sun)) -pi/2; 
eta_angle = cos(a_sun);
eta_angle(isnan(eta_angle)) = cosd(4);

% degradation
eta_degradation = (1-D).^t;

%total efficency
eta = eta_space.*eta_angle.*eta_degradation;

eta = eta/max(eta);
T_max_modulated = T_max*eta;



end






