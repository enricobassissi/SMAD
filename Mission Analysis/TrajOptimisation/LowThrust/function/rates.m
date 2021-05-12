function dydt = rates(t,f,id)
% ------------------------
%{
This function calculates the acceleration vector using Equation 2.22.
t - time
f - column vector containing the position vector and the
velocity vector at time t
x, y, z - components of the position vector r
r - the magnitude of the position vector
vx, vy, vz - components of the velocity vector v
ax, ay, az - components of the acceleration vector a
dydt - column vector containing the velocity and acceleration
components
id - name of celestial object around which the s/c is orbiting, to be used
in the switch loop to find the right Gravitational constant


CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John

%}
% ------------------------

 switch id
          case 'sun'
              mu=1.327124400180000e+11;
          case 'earth'
              mu =  3.98600433e+5;
          case 'mercury'
              mu=2.203208e4;
          case 'mars'
              mu=4.2828314E+4;
          case 'neptune'
              mu=6.83653406400E+06;
          case 'moon'
              mu = astroConstants(20);
          otherwise
            disp('Unknown gravitational constant')
 end

x = f(1);
y = f(2);
z = f(3);
vx = f(4);
vy = f(5);
vz = f(6);
r = norm([x y z]);
ax = -mu*x/r^3;
ay = -mu*y/r^3;
az = -mu*z/r^3;
dydt = [vx vy vz ax ay az]';
end 