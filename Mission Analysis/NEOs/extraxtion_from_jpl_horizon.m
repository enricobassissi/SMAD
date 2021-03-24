function [r, v, sim] = extraxtion_from_jpl_horizon(table_name)

%{
INPUT:
    table_name: string

OUTPUT:
    r: position vector, km
    v: velocity vector, km/s
    sim: structure with infos extracted from the jpl nasa horizon simulation

FUNCTIONS NEEDED:
    sv_from_coe.m
    astroConstants.m
%}

% Data extraction from table of EROS ephemeris
data_table = readtable(table_name);
data_array = str2double(table2array(data_table(:,5:end)));

% Astronomical constants
AU = astroConstants(2);
muSun = astroConstants(4);

% Indexes corresponding to the correct order of the keplerian elements
% inside the data: [a,e,i,OM,om,th]
idx_coe_in_data_table = [10,1,3,4,5,9]; % [AU,-,deg,deg,deg,deg]
coe = [data_array(:,idx_coe_in_data_table(1))*AU,...
            data_array(:,idx_coe_in_data_table(2)),...
            deg2rad(data_array(:,idx_coe_in_data_table(3))),...
            deg2rad(data_array(:,idx_coe_in_data_table(4))),...
            deg2rad(data_array(:,idx_coe_in_data_table(5))),...
            deg2rad(data_array(:,idx_coe_in_data_table(6)))]; % [km,-,rad,rad,rad,rad]

% Generation of corresponding state vectors
r = zeros(length(coe),3);
v = zeros(length(coe),3);
for i = 1:length(coe)
    [r(i,:), v(i,:)] = sv_from_coe(coe(i,:), muSun);
end

% Date of the ephemeris analysed
sim.date_vect = datevec(table2array(data_table(:, 3)));
sim.date_start = datevec(table2array(data_table(1, 3)));
sim.date_end = datevec(table2array(data_table(end, 3)));
sim.mjd2000_start = date2mjd2000(sim.date_start);
sim.mjd2000_end = date2mjd2000(sim.date_end);

end