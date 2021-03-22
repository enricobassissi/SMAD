function [r, v] = sv_from_coe(coe,mu)

%{
This function computes the state vector (r,v) from the classical orbital elements (coe).
 INPUT: 1. mu  : gravitational parameter (km^3;s^2)
        2. coe : orbital elements [a, e, incl, OM, om, f]

 OUTPUT: 1. r  : position vector in the geocentric equatorial frame (km)
         2. v  : velocity vector in the geocentric equatorial frame (km/s)
  
 FUNCTIONS REQUIRED: -

 CONTRIBUTORS:  Bassissi Enrico
                Colombo Alessandro
                De Luca Maria Alessandra

 QUANTITIES INVOLVED: 1. h : angular momentum (km^2/s)
                      2. e : eccentricity
                      3. OM : right ascension of the ascending node (rad)
                      4. incl : inclination of the orbit (rad)
                      5. om : argument of perigee (rad)
                      6. f : true anomaly (rad)
                      7. R3_w : Rotation matrix about the z-axis through the angle w
                      8. R1_i : Rotation matrix about the x-axis through the angle i
                      9. R3_W : Rotation matrix about the z-axis through the angle RA
                      10. Q_pX : Matrix of the transformation from perifocal to geocentric equatorial frame
                      11. rp : position vector in the perifocal frame (km)
                      12. vp : velocity vector in the perifocal frame (km/s)
%}

% Inizialization
a = coe(1);
e = coe(2);
incl = coe(3);
OM = coe(4);
om = coe(5);
f = coe(6);

h = sqrt(mu*a*(1-e^2));

rp = (h^2/mu) * (1/(1 + e*cos(f))) * (cos(f)*[1;0;0] + sin(f)*[0;1;0]);
vp = (mu/h) * (-sin(f)*[1;0;0] + (e + cos(f))*[0;1;0]);

% First rotation of OM (right ascension of the ascending node (rad))
R3_OM = [ cos(OM) sin(OM) 0
        -sin(OM) cos(OM) 0
            0 0 1];
        
% Second rotation of incl ( inclination of the orbit (rad))
R1_i = [1 0 0
        0 cos(incl) sin(incl)
        0 -sin(incl) cos(incl)];
    
% Third rotation of om (argument of perigee (rad))
R3_om = [ cos(om) sin(om) 0
    -sin(om) cos(om) 0
    0 0 1];

% Final position, rotation combined
Q_pX = (R3_om*R1_i*R3_OM)';

% Final vectors position and velocity (r and v are column vectors):
r = Q_pX*rp;
v = Q_pX*vp;

end