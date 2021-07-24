function coasting = time_and_distance_at_asteroid(ast_name, arr_mjd2000, dep_mjd2000, data, sim, n)

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
coasting.v_ast = v_ast;
coasting.rEA = rEA;
coasting.vEA = vEA;

% ---- frame centered in the earth ------ %
% ---- relative positions asteroid wrt earth ----- %
coasting.rel_pos_ast_earth =  r_ast/sim.DU - rEA/sim.DU;
coasting.norm_rel_pos_ast_earth = vecnorm(coasting.rel_pos_ast_earth,2,2);

coasting.rel_vel_ast_earth =  v_ast/sim.DU*sim.TU - vEA/sim.DU*sim.TU;
coasting.norm_rel_vel = vecnorm(coasting.rel_vel_ast_earth,2,2);


% ---- frame centered in the asteroid ------ %
% ---- relative positions earth wrt asteroid ----- %
% -- but still not rotating with the asteroid
coasting.rel_pos_earth_ast =  rEA/sim.DU - r_ast/sim.DU;
coasting.norm_rel_pos_earth_ast = vecnorm(coasting.rel_pos_earth_ast,2,2);

% ---- frame centered in the asteroid ------ %
% ---- spherical coordinates of the earth as seen from the ast ----- %
% ---- frame moving with the asteroid ---- %
% versor of asteroid position, point by point, positive outward
r_ast_norm = vecnorm(r_ast,2,2);
r_ast_vers = r_ast./r_ast_norm;
% versor of velocity , clockwise
v_ast_norm = vecnorm(v_ast,2,2);
v_ast_vers = v_ast./v_ast_norm;
% versor normal to the plane of the ast orbit
h_ast_vers = cross(r_ast_vers,v_ast_vers);
% tg versor, normal to r_vers and h_vers
tg_ast_vers = cross(r_ast_vers,h_ast_vers);

% projection of rel position of the earth on the frame changing with the asteroid
x_EA_frame_asteroid = zeros(n,3);
y_EA_frame_asteroid = zeros(n,3);
z_EA_frame_asteroid = zeros(n,3);
for i= 1:n
    dot_prod_x = dot(coasting.rel_pos_earth_ast(i,:),r_ast_vers(i,:));
    dot_prod_y = dot(coasting.rel_pos_earth_ast(i,:),tg_ast_vers(i,:));
    dot_prod_z = dot(coasting.rel_pos_earth_ast(i,:),h_ast_vers(i,:));
    x_EA_frame_asteroid(i,:) = dot_prod_x*r_ast_vers(i,:);
    y_EA_frame_asteroid(i,:) = dot_prod_y*tg_ast_vers(i,:);
    z_EA_frame_asteroid(i,:) = dot_prod_z*h_ast_vers(i,:);
end

coasting.position_EA_frame_asteroid = x_EA_frame_asteroid+y_EA_frame_asteroid+z_EA_frame_asteroid;

% position_EA_frame_asteroid_norm = vecnorm(position_EA_frame_asteroid,2,2);
% TOL = position_EA_frame_asteroid_norm - coasting.norm_rel_pos_earth_ast;
% max(abs(TOL)) % AU


end