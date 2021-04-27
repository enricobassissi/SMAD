function array_horizons = read_web_form_of_JPL_Horizons(txt_name_string_dot_txt)
% JDTDB    Julian Day Number, Barycentric Dynamical Time
% EC (5)     Eccentricity, e 
% QR     Periapsis distance, q (au)
% IN (7)     Inclination w.r.t X-Y plane, i (degrees)
% OM (8)    Longitude of Ascending Node, OMEGA, (degrees)
% W  (9)    Argument of Perifocus, w (degrees)
% Tp     Time of periapsis (Julian Day Number)
% N      Mean motion, n (degrees/day)
% MA     Mean anomaly, M (degrees)
% TA (13)    True anomaly, nu (degrees)
% A  (14)    Semi-major axis, a (au)
% AD     Apoapsis distance (au)
% PR     Sidereal orbit period (day)
% 

temp = readtable(txt_name_string_dot_txt);
array_horizons = str2double(table2array([temp(:,14),temp(:,5),temp(:,7),temp(:,8),temp(:,9),temp(:,13)]));

end