function [Gd] = solar_irradiance_cooler(d, Gref, dref)
%	soloar_irradiance - Compute the solar irradiance at a certain distance
%	from the sun CONSIDERING (very roughly) THE INCREASE IN EFFICIENCY OF THE PANLES DUE
%	TO A LOWER TEMPERATURE AT HIGHER DISTANCES.
%
% PROTOTYPE:
%  [Gd] = solar_irradiance_cooler(d, Gref, dref)
%
% DESCRIPTION:
%   Compute the solar irradiance "Gd" at distance"d", knowing a reference 
%   solar irradiacne "Gref" at distance "dref", using the inverse square law. 
%
%  INPUT :
%   d:      [1] Distance from the sun center, where to evaluate the irradiance
%   Gref:   [1] Solar irradiance at reference distance, usually solar
%               constant = 1361 W/m2 ECSS-E-ST-10-04C-Rev.1 (15June2020) 
%   dref:   [1] Reference distance, usually = 1 AU
%
%  OUTPUT:
%   Gd:     [1] Solar irradiance at selected distance
%
%  FUNCTIONS CALLED:
%
% AUTHOR:
%   Marco Elia
%
% PREVIOUS VERSION:
%   \\
%
% CHANGELOG:
%   20/05/2021, Marco Elia
%
% -------------------------------------------------------------------------
% REFERENCES:
%   SPACECRAFT SYSTEMS ENGINEERING Fourth Edition Edited by Peter Fortescue
%   pp 341
Gd = Gref*(dref./d).^1.5;
end

