function [coe] = coe_from_sv(R,V,mu)
% -----------------------------------
%{
% This function computes the classical orbital elements (coe)...
% from the state vector (R,V)
% coe - vector of orbital elements [a, e, incl, OM, om, f, h]...
% as dimension [km, /, rad, rad, rad, rad, km^2/s]
mu - gravitational parameter (km^3/s^2)
R - position vector in the geocentric equatorial frame (km)
V - velocity vector in the geocentric equatorial frame (km)
r, v - the magnitudes of R and V
vr - radial velocity component (km/s)
H - the angular momentum vector (km^2/s)
h - the magnitude of H (km^2/s)
incl - inclination of the orbit (rad)
N - the node line vector (km^2/s)
n - the magnitude of N
cp - cross product of N and R
OM - right ascension of the ascending node (rad)
E - eccentricity vector
e - eccentricity (magnitude of E)
eps - a small number below which the eccentricity is considered
to be zero
om - argument of perigee (rad)
f - true anomaly (rad)
a - semimajor axis (km)
User M-functions required: None

CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John

%}
% ---------------------------------------------
% 00. intro parameter, tolerance
eps = 1.e-10;
% 0. magnitudes of R and V
r = norm(R);
v = norm(V);
% 1. radial velocity component (km/s)
vr = dot(R,V)/r;
% 2. angular momentum vector (km^2/s) and magnitude
H = cross(R,V);
h = norm(H);
% 3. inclination of the orbit (rad)
incl = acos(H(3)/h);
% 4. the node line vector (km^2/s) and magnitude
N = cross([0 0 1],H);
n = norm(N);
% 5. right ascension of the ascending node (rad)
if n ~= 0
OM = acos(N(1)/n);
if N(2) < 0
OM = 2*pi - OM;
end
else
OM = 0;
end
% 6. eccentricity vector
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);
% 7. argument of perigee (rad)  (incorporating the case e = 0):
if n ~= 0
    if e > eps
        om = acos(dot(N,E)/n/e);
        if E(3) < 0
            om = 2*pi - om;
        end
    else
        om = 0;
    end
else
    om = 0;
end
% 8. true anomaly (rad) (incorporating the case e = 0):
if e > eps
    f = acos(dot(E,R)/e/r);
    if vr < 0
        f = 2*pi - f;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        f = acos(dot(N,R)/n/r);
    else
        f = 2*pi - acos(dot(N,R)/n/r);
    end
end
% 9. semimajor-axis (km) (observe that a < 0 for a hyperbola):
a = h^2/mu/(1 - e^2);
coe = [a, e, incl, OM, om, f, h];
end