function [TOF, T] = time_law(f1,f2,a,e,mu)
%{
VALID ONLY FOR 1 REVOLUTION EXTRA, K=1
a [km]
f [rad]
e []
mu [km^3/s^2]

TOF [s]
T_half [s]
%}

%{
example for validation:
-- inputs --
theta1 = 30°;
theta2 = 90°;
a = 7500; km
e = 0.1;
mu_E = 398600;
[dT] = time_law(deg2rad(30),deg2rad(90),7500,0.1,mu_E)
-- procedures results --
t0 = 0;
E1 = 0.47557 radians
E2 = 1.47063 radians
M1 = 0.42978 radians
M2 = 1.37113 radians
n = 0.00097202 rad/s
t = 0 + (1.37113 - 0.42978) / 0.00097202
-- final results! --
t = 968.4 s

note that if f2 < f1, the result will be negative!
the result is robust even if going over 360°
%}

t0 = 0;

% cos E = (e + cos nu0) / (1 + e cos )
E1 = acos((e + cos(f1)) / (1 + e*cos(f1)));
E2 = acos((e + cos(f2)) / (1 + e*cos(f2)));

% Some care is required if the spacecraft passes k times through periapsis, as 2kπ must be added to E2.
if f1 > pi
    E1 = 2*pi - E1;
end
    
if f2 > pi && f2 <= 2*pi
    E2 = 2*pi - E2;
elseif f2 > 2*pi && f2 <= 3*pi
    E2 = E2 + 2*pi;
elseif f2 > 3*pi
    E2 = 4*pi - E2;
end

% M = E - e × sin E
M1 = E1 - e*sin(E1);
M2 = E2 - e*sin(E2);

% n = sqrt[ GM / a3 ] || n = sqrt( mu / a^3 );
n = sqrt( mu / a^3 );

% M - M0 = n × (t - t0)
TOF = t0 + (M2 - M1) / n;

T = 2*pi/n; % period of that orbit

end