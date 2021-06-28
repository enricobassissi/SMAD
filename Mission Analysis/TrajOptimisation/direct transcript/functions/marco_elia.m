%% --- stuff per marco elia
function [EPS] = marco_elia(SC, r_encounter, data, sim, href, colors)

EPS.time = SC.timead;

% ----- cartesian 3D trajectory
r = SC.HS.X(:,1);
theta = SC.HS.X(:,2);
z = SC.HS.X(:,3);
vr = SC.HS.X(:,4);
vtheta = SC.HS.X(:,1).*SC.HS.X(:,5);
vz = SC.HS.X(:,6);

r_transf_orbit  = [r.*cos(theta), r.*sin(theta), z];
EPS.R_cartesian = rotate_local2ecplitic(r_encounter,r_transf_orbit,data.n_int,href);

% --- local cartesian
[x,y,z] = pol2cart(theta,r,z);
in_plane_position = [x,y,zeros(length(x),1)];
in_plane_position_norm = vecnorm(in_plane_position,2,2);
in_plane_position_vers = in_plane_position./in_plane_position_norm;

z_vers = [0,0,1].*ones(length(x),3);
tg_vers = -cross(in_plane_position_vers,z_vers);

% ----- cartesian 3D thrust
% beta = elevation, alpha = thrust angle in plane
% x -> T t , y = -T rad, z -> T outofplane
[Tx, Ty, Tz] = sph2cart(SC.HS.alpha,SC.HS.beta,SC.HS.T);
% cartesiano centrato movente con la spacecraft
Ttg = Tx;
Tradial = -Ty;
Tout = Tz;

T_cart_radial = Tradial.*in_plane_position_vers; 
T_cart_tg = Ttg.*tg_vers;
T_cart_out = Tout.*z_vers;

T_cart_local = T_cart_radial+T_cart_tg+T_cart_out;

% check -- if low it's the cross product
% find(SC.HS.T == vecnorm(T_cart_local,2,2));
% min(SC.HS.T - vecnorm(T_cart_local,2,2));
% max(SC.HS.T - vecnorm(T_cart_local,2,2));

% -- T in cartesiano
EPS.T_transf_orbit = rotate_local2ecplitic(r_encounter,T_cart_local,data.n_int,href);

% --- plot orbit 3D
figure('Name','Orbit Plot')
plot3(EPS.R_cartesian(:,1)./sim.DU,EPS.R_cartesian(:,2)./sim.DU,EPS.R_cartesian(:,3)./sim.DU,...
    'Color',colors(1,:));
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');

% -- plot thrust quiver
figure('Name','Thrust Plot')
plot3(EPS.R_cartesian(:,1)./sim.DU,EPS.R_cartesian(:,2)./sim.DU,EPS.R_cartesian(:,3)./sim.DU,...
    'Color',colors(1,:));
hold on
quiver3(EPS.R_cartesian(:,1)./sim.DU,EPS.R_cartesian(:,2)./sim.DU,EPS.R_cartesian(:,3)./sim.DU,EPS.T_transf_orbit(:,1),...
    EPS.T_transf_orbit(:,2),EPS.T_transf_orbit(:,3),0.5)
axis equal; grid on;
xlabel('x [AU]'); ylabel('y [AU]'); ylabel('y [AU]'); 
title('In-plane + out-of-plane Thrust')
% Sun
plot3(0,0,0,'o','Color',colors(4,:),'DisplayName','Sun')
% legend('show')
view(2)

end