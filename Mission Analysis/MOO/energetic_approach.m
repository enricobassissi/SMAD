clear; close; clc;
% Import module of Python
module = py.importlib.import_module('neo_api_function');
AU = astroConstants(2); %km
muSun = astroConstants(4); %km^3/s^2

%% sc data
% dawn, hyabusa 2,deep space 1
m = [740, 490, 373]; %kg 
T = [90, 28*3, 90]; %mN
Isp = [3100, 3000,3280]; %s
P = [1e4, 1250, 2300]; %W

m_our = 100; %kg
T_our = interp1(m, T, m_our, 'spline'); %interpolate data to obtain conductivity

figure()
plot(m,T)
hold on
plot(m_our,T_our,'*')

%% energetic approach
m_dry = 373; % kg, dry
T_our = 90e-3; %N
acc = T_our/m_dry; %m/s^2
% Earth at dep
mjd0 = date2mjd2000([1998 10 24 0 0 0]);
[kep_EA,ksun] = uplanet(mjd0, 3);
[r1, ~] = sv_from_coe(kep_EA,ksun);
r_in = norm(r1); %km
% Ast at arr
PointOfView = 'Sun';
epoch_start = '1999-07-28';
epoch_stop = '1999-07-29';
step = '1d';
type_elements = 'Vectors';
ast_name = '1992KD';
% data extraction section
py_data = py.neo_api_function.get_horizons_ephemerides(py.str(ast_name),py.str(PointOfView),...
                  py.str(epoch_start),py.str(epoch_stop),py.str(step),py.str(type_elements));
horizons_data = double(py_data);

r_fin = norm(horizons_data(end, 1:3))*AU; %km

TOF = 1/acc*((-sqrt(muSun/r_fin) + sqrt(muSun/r_in))*1e3); %s 
dV = acc*TOF; % m/s
g0 = 9.81; %m/s^2
Isp = 3280; %s
m_prop = m_dry*(exp(dV/(g0*Isp)) - 1); %kg

% print table values
Quantity = {'dV tot [km/s]';'TOF tot [d]';'m prop [kg]'};
table(Quantity,[dV*1e-3;TOF/(3600*24);m_prop])

