function [r, v] = sv_from_coe(coe,mu)
% ---------------------------------
%{
This function computes the state vector (r,v) from the
classical orbital elements (coe).
mu - gravitational parameter (km^3;s^2)
coe - orbital elements [a, e, incl, OM, om, f, h]
where
h = angular momentum (km^2/s)
e = eccentricity
OM = right ascension of the ascending node (rad)
incl = inclination of the orbit (rad)
om = argument of perigee (rad)
f = true anomaly (rad)
R3_w - Rotation matrix about the z-axis through the angle w
R1_i - Rotation matrix about the x-axis through the angle i
R3_W - Rotation matrix about the z-axis through the angle RA
Q_pX - Matrix of the transformation from perifocal to geocentric
equatorial frame
rp - position vector in the perifocal frame (km)
vp - velocity vector in the perifocal frame (km/s)
r - position vector in the geocentric equatorial frame (km)
v - velocity vector in the geocentric equatorial frame (km/s)
User M-functions required: none


CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John

%}
% ----------------------------------------------
a = coe(1);
e = coe(2);
incl = coe(3);
OM = coe(4);
om = coe(5);
f = coe(6);
p=a*(1-e^2);
h=sqrt(p*mu);



% note that rp and vp are column vectors:
rp = (h^2/mu) * (1/(1 + e*cos(f))) * (cos(f)*[1;0;0] + sin(f)*[0;1;0]);
vp = (mu/h) * (-sin(f)*[1;0;0] + (e + cos(f))*[0;1;0]);
% first rotation of OM (right ascension of the ascending node (rad))
R3_OM = [ cos(OM) sin(OM) 0
        -sin(OM) cos(OM) 0
            0 0 1];
% second rotation of incl ( inclination of the orbit (rad))
R1_i = [1 0 0
        0 cos(incl) sin(incl)
        0 -sin(incl) cos(incl)];
% third rotation of om (argument of perigee (rad))
R3_om = [ cos(om) sin(om) 0
    -sin(om) cos(om) 0
    0 0 1];
% final position, rotation combined
Q_pX = (R3_om*R1_i*R3_OM)';
% final vectors position and velocity (r and v are column vectors):
r = Q_pX*rp;
v = Q_pX*vp;

end