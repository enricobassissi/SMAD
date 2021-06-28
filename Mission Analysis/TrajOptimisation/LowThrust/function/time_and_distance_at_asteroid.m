function coasting = time_and_distance_at_asteroid(ast_name, arr_mjd2000, dep_mjd2000, data, sim)

n = sim.n_sol;

time_vect = linspace(arr_mjd2000,dep_mjd2000,n);

rEA = zeros(n,3); vEA = zeros(n,3);
r_ast = zeros(n,3); v_ast = zeros(n,3);
kep_EA = zeros(n,6); kep_ast = zeros(n,6);
for i = 1:n
    [kep_EA(i,:), ksun] = uplanet(time_vect(i), 3); % Earth
    [rEA(i,1:3), vEA(i,1:3)] = sv_from_coe(kep_EA(i,:),ksun); %km, km/s
    
    [kep_ast(i,:)] = uNEO3(time_vect(i),ast_name,data);
    [r_ast(i,1:3),v_ast(i,1:3)] = sv_from_coe(kep_ast(i,:),ksun); % km, km/s
end

coasting.time = time_vect';

coasting.coe_ast = kep_ast;
coasting.coe_EA = kep_EA;

% ---- frame centered in the sun ------ %
% ---- heliocentric absolute positions ----- %
coasting.r_ast = r_ast;
coasting.rEA = rEA;

% ---- frame centered in the earth ------ %
% ---- relative positions asteroid wrt earth ----- %
coasting.rel_pos_ast_earth =  r_ast/sim.DU - rEA/sim.DU;
coasting.norm_rel_pos_ast_earth = vecnorm(coasting.rel_pos_ast_earth,2,2);

coasting.rel_vel_ast_earth =  v_ast/sim.DU*sim.TU - vEA/sim.DU*sim.TU;
coasting.norm_rel_vel = vecnorm(coasting.rel_vel_ast_earth,2,2);


% ---- frame centered in the asteroid ------ %
% ---- relative positions earth wrt asteroid ----- %
coasting.rel_pos_earth_ast =  rEA/sim.DU - r_ast/sim.DU;
coasting.norm_rel_pos_earth_ast = vecnorm(coasting.rel_pos_earth_ast,2,2);

[coasting.az,coasting.el,coasting.r] = cart2sph(coasting.rel_pos_earth_ast);


end