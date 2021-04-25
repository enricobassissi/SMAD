%% ----- Analytic Approach ----- %%
%% Introduction
clear; close; clc;
% Import module of Python
module = py.importlib.import_module('neo_api_function');

addpath time
addpath function
AU = astroConstants(2); %km
muSun = astroConstants(4); %km^3/s^2

%% Calculations
m_dry = 373; % kg, dry
m_wet = 373+80; %kg
T_our = 90e-3; %N
acc = T_our/m_wet; %m/s^2
% Earth at dep
mjd0 = date2mjd2000([1998 10 24 0 0 0]);
[kep_EA,ksun] = uplanet(mjd0, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);
r_in = norm(r1); %km
v_in = norm(v1); %km/s
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
%%
r_fin = norm([-6.19009882e-01 -1.17474220e+00  7.10312086e-04].*AU);

% analytic approach stuff
TOF = v_in*1e3/acc*(1 - (2*acc*r_in*1e3/(v_in*1e3)^2)^1/4); %s
TOF_d = TOF/(3600*24)
TOF_true = date2mjd2000([1999 07 29 0 0 0]) - mjd0;
%%
dV = acc*TOF; % m/s
dV_kms = dV*1e-3
g0 = 9.81; %m/s^2
Isp = 3280; %s
m_prop = m_dry*(exp(dV/(g0*Isp)) - 1); %kg

%% print table values
Quantity = {'dV tot [km/s]';'TOF tot [d]';'m prop [kg]'};
table(Quantity,[dV*1e-3;TOF/(3600*24);m_prop])


%% TERRA MARTE ESEMPIO PRINETTO
clear; clc

addpath time
addpath function
AU = astroConstants(2); %km
muSun = astroConstants(4); %km^3/s^2

%% Calculations
m_dry = 1000; % kg, dry
% m_wet = 373+80; %kg
T_our = 220e-3; %N
acc = T_our/m_dry; %m/s^2
% Earth at dep
mjd0 = date2mjd2000([2030 05 01 0 0 0]);
[kep_EA,ksun] = uplanet(mjd0, 3);
[r1, v1] = sv_from_coe(kep_EA,ksun);
r_in = norm(r1)*1e3; %m
v_in = norm(v1)*1e3; %m/s
% Earth at dep
mjdf = date2mjd2000([2032 05 01 0 0 0]);
[kep_MA,ksun] = uplanet(mjdf, 4);
[r2, v2] = sv_from_coe(kep_MA,ksun);
r_fin = norm(r2)*1e3; %m
v_fin = norm(v2)*1e3; %m/s

% analytic approach stuff
% TOF = v_in*1e3/acc*(1 - (2*acc*r_in*1e3/(v_in*1e3)^2)^1/4); %s
TOF = (1-r_in/r_fin)/(acc*sqrt(r_in/(ksun*1e3)))
TOF_d = TOF/(3600*24)
% TOF_true = date2mjd2000([1999 07 29 0 0 0]) - mjd0;

dV = acc*TOF; % m/s
dV_kms = dV*1e-3
g0 = 9.81; %m/s^2
Isp = 3000; %s
m_prop = m_dry*(exp(dV/(g0*Isp)) - 1); %kg
