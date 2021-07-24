function [v_cart_local] = from_dt_sol_to_cart(X)

vr = X(:,4); % 4 vr, 5 theta_dot, 6 vz
vth = X(:,1).*X(:,5);
vz = X(:,6);

% --- local cartesian
r = X(:,1);
theta = X(:,2);
z = X(:,3);
[x,y,z] = pol2cart(theta,r,z);
in_plane_position = [x,y,zeros(length(x),1)];
in_plane_position_norm = vecnorm(in_plane_position,2,2);
in_plane_position_vers = in_plane_position./in_plane_position_norm;

z_vers = [0,0,1].*ones(length(x),3);
tg_vers = -cross(in_plane_position_vers,z_vers);

v_cart_radial = vr.*in_plane_position_vers; 
v_cart_tg = vth.*tg_vers;
v_cart_out = vz.*z_vers;

v_cart_local = v_cart_radial+v_cart_tg+v_cart_out;

end