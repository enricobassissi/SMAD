function [T_available] = available_thrust(t, RR, TT, data)
%	available_power_thrust - Compute the thrust available during the mission
%
% PROTOTYPE:
%  [P_propulsion] = available_power_for_propulsion(t, RR, TT, data)
%
% DESCRIPTION:
%   Compute the power available for the
%	propulsion system during the mission
%
%  INPUT :
%   t:      [Nx1] Array of times [years]
%   RR:     [Nx3] Array of vectors position in heliocentric cartesian frame [AU]
%   TT:     [Nx3] Array of vectors thrust (needed just for the orientation)
%   data:   Parameters referrred to thruster solar panels etc...
%         P_sa:                         [1] Power generated by the solar array [W]
%         eta_one_string_failure:       [1] Efficency due to one string failure [-]
%         eta_temperature:              [1] Rfficiency due to temperature [-]
%         D_cell_radiation:             [1] Degradation of cells efficency due to radiation per year [%/year]
%         D_cell_other:                 [1] Degradation of panels efficency due to other effects per year [%/year]
%         P_other_subsystem_margined:	[1] Power margined used by othger subsystems [W]
%         X_d:                          [1] Efficiency due to DET power transfer [-]
%         component_margin:             [1] Power margin on the component [-] 1.05=+5%
%         a_fitted:                     [1] Coeff of linear fitting of Power thrust curve P = a*T+b
%         b_fitted:                     [1] Coeff of linear fitting of Power thrust curve P = a*T+b
%         T_limit:                      [1] Maximum thrust that can be generated by the thruster [mN]
%
%  OUTPUT:
%   T_available:                       [Nx1] Available thrust [mN]
%
%  FUNCTIONS CALLED:
%   solar_irradiance_cooler
%   available_power_for_propulsion
%
% AUTHOR:
%   Marco Elia
%
% PREVIOUS VERSION:
%   \\
%
% CHANGELOG:
%   07/06/2021, Marco Elia
%
% -------------------------------------------------------------------------
% REFERENCES:

P_propulsion = available_power_for_propulsion(t, RR, TT, data);

T_available = (P_propulsion-data.b_fitted)/data.a_fitted;% computed from fitting below 
T_available(T_available > data.T_limit) = data.T_limit; % capping thrust
end






