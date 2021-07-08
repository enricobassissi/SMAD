%% MOVIE PLOT OF THE WHOLE MISSION SC2
%{
and the full orbit of inner solar system in that same 
time of the mission
%}
%% Default options
clear; close all; clc;

set(0, 'DefaultTextFontSize', 20) % modify it if too small
set(0, 'DefaultAxesFontSize', 20) % modify it if too small
set(0, 'DefaultLegendFontSize', 20) % modify it if too small
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.8)

colors = [0    50   71;... % (1) DEEP SPACE
          207  29   57;... % (2) EXCITE RED +1
          0    142  122;... % (3) PURE TEAL +1
          251  171  24;... % (4) ENLIGHT YELLOW
          244  121  32;... % (5) ENLIGHT YELLOW +1
          150  1    54;... % (6) EXCITE RED +2
          167  85   52;... % (7) ENLIGHT YELLOW +2
          0    97   158;... % (8) TRUSTY AZURE +1
          30   51   120;... % (9) TRUSTY AZURE +2
          0    103  98;... % (10) PURE TEAL +2
          51   94   111;... % (11) DEEP SPACE -1
          0    0    0]./255; % (12) BLACK

%% work environment setup
str_path=split(pwd, 'PostAnalysis\Video_drawnow');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
dt_path=string(str_path(1))+'TrajOptimisation\direct transcript\follower';
addpath(genpath(dt_path));
fun_path=string(str_path(1))+'TrajOptimisation\direct transcript\functions';
addpath(genpath(fun_path));

%% load and textures
load('Mission_160kg_dry.mat')
load('MA5_160kg_dry.mat')
load('data_elements_matrix_44_63_2SC.mat')

%% --------------------------------------------------------------------- %%
%% --------------------------------- SC2 ------------------------------- %%
%% --------------------------------------------------------------------- %%
%% the time vector
% MA.SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
%     SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000
time_vect = MA.SC2.mjd2000_dep+MA.SC2.uniform.time/86400;
err_1st_tof = time_vect(101)-time_vect(100)
new_tofa = time_vect(100)-time_vect(1)
new_cta = SC1.CT1+time_vect(101)-time_vect(100)
err_2nd_tof = time_vect(301)-time_vect(300)
new_tofb = time_vect(300)-time_vect(201)

% time_vect2 = MA.SC2.mjd2000_dep+MA.SC2.uniform.time/86400;
% time_vect2(101)-time_vect2(100)
% time_vect2(301)-time_vect2(300)
% -- length
n = length(time_vect);

%% sc trajectory
% --- new encounters
MJD01_dim = SC2.mjd2000_dep;
[kep_EA,~] = uplanet(MJD01_dim, 3);
[r_EA, v_EA] = sv_from_coe(kep_EA,sim.mu_dim);
r_EA = r_EA/sim.DU;
v_EA = v_EA/sim.DU*sim.TU;
new_r_encounter.EA = r_EA;
new_v_encounter.EA = v_EA;

MJDAa_dim = (SC2.mjd2000_dep+new_tofa);
[kep_ast_Aa] = uNEO3(MJDAa_dim,SC2.asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rAa, vAa] = sv_from_coe(kep_ast_Aa,sim.mu_dim); % km, km/s
rAa = rAa/sim.DU;
vAa = vAa/sim.DU*sim.TU;
new_r_encounter.astAa = rAa;
new_v_encounter.astAa = vAa;

MJDDa_dim = (SC2.mjd2000_dep+new_tofa+new_cta);
[kep_ast_Da] = uNEO3(MJDDa_dim,SC2.asteroid_a,data); % [km,-,rad,rad,rad,wrapped rad]
[rDa, vDa] = sv_from_coe(kep_ast_Da,sim.mu_dim); % km, km/s
rDa = rDa/sim.DU;
vDa = vDa/sim.DU*sim.TU;
new_r_encounter.astDa = rDa;
new_v_encounter.astDa = vDa;

MJDAb_dim = (SC2.mjd2000_dep+new_tofa+new_cta+new_tofb);
[kep_ast_Ab] = uNEO3(MJDAb_dim,SC2.asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rAb, vAb] = sv_from_coe(kep_ast_Ab,sim.mu_dim); % km, km/s
rAb = rAb/sim.DU;
vAb = vAb/sim.DU*sim.TU;
new_r_encounter.astAb = rAb;
new_v_encounter.astAb = vAb;

MJDDb_dim = (SC2.mjd2000_dep+new_tofa+new_cta+new_tofb+new_cta);
[kep_ast_Db] = uNEO3(MJDDb_dim,SC2.asteroid_b,data); % [km,-,rad,rad,rad,wrapped rad]
[rDb, vDb] = sv_from_coe(kep_ast_Db,sim.mu_dim); % km, km/s
rDb = rDb/sim.DU;
vDb = vDb/sim.DU*sim.TU;
new_r_encounter.astDb = rDb;
new_v_encounter.astDb = vDb;

%% DT - SC 2 - LEG a
data.Tmax = 0.015; % max thrust available for EGI
data.Pmax = 450; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2750/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 2;
data.angle_inplane_panels_max = 10; % degrees
data.n_int = 50;
data.muS = sim.mu_dim;   
data.MU = sol.M_start_SC2_lega; % mass adimensionalisation on wet mass
SC2.asteroid_a = sol.asteroid_a;

v_launcher = sol.v_inf_magn/sim.DU*sim.TU*[cos(sol.el)*cos(sol.az); cos(sol.el)*sin(sol.az); sin(sol.el)];
v_dep = new_v_encounter.EA + v_launcher;  %if parabolic escape (v_extra = 0)
SC2.lega.v_in = v_dep;
SC2.lega.v_end = new_v_encounter.astAa;
SC2.lega.r_in = new_r_encounter.EA;
SC2.lega.r_end = new_r_encounter.astAa;
SC2.lega.N_rev = sol.Nrev(3);
SC2.lega.TOF = new_tofa*86400/sim.TU;
SC2.lega.Href = sol.Href(:,3);

DT_input = SC2.lega;
tic
[SC2.lega.HS, SC2.lega.RES, SC2.lega.timead, SC2.lega.Xad, SC2.lega.Xpropad, ...
    SC2.lega.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC2_lega = toc/60; % minutes

SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.lega.EPS = marco_elia(SC2.lega,new_r_encounter.EA,data,sim,SC2.lega.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC2.lega] = interpoliamo_porc(SC2.lega,data.n_int);
SC2.lega.v_in = v_dep;
SC2.lega.v_end = new_v_encounter.astAa;
SC2.lega.r_in = new_r_encounter.EA;
SC2.lega.r_end = new_r_encounter.astAa;
SC2.lega.N_rev = sol.Nrev(3);
SC2.lega.TOF = new_tofa*86400/sim.TU;
SC2.lega.Href = sol.Href(:,3);
SC2.mjd2000_dep = sol.departure_mjd2000;

SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.lega.EPS = marco_elia(SC2.lega,new_r_encounter.EA,data,sim,SC2.lega.Href,colors);
% -- new mass after DT calculations
SC2.lega.mass_start = SC2.lega.HS.X(1,7);
SC2.lega.mass_end = SC2.lega.HS.X(end,7);
SC2.lega.mass_depleted = SC2.lega.mass_start - SC2.lega.mass_end;
SC2.lega.mass_fraction = SC2.lega.mass_depleted/SC2.lega.mass_start;

% -- coasting between the two legs
[SC2.coasting.lega] = relative_position_coasting_stuff(SC2.asteroid_a,SC2.mjd2000_dep,...
    new_tofa,new_cta,data.n_int,EphData,sim,colors);

%% DT - SC 2 - LEG b
data.Tmax = 0.02; % max thrust available for EGI
data.Pmax = 550; % max power related to 15 mN thrust on qitetiq t5
data.Isp = 2900/sim.TU; % for 15 mN the Isp of Qinetiq T5 is 2750s
data.ThrustModulationFlag = 1;
data.angle_inplane_panels_max = 20; % degrees
data.n_int = 20;
% data.MU = sol.M_start_SC2_legb; % mass adimensionalisation on wet mass
data.MU = SC2.lega.mass_end - sim.M_pods; % mass adimensionalisation on wet mass
SC2.asteroid_b = sol.asteroid_b;

SC2.legb.v_in = new_v_encounter.astDa;
SC2.legb.v_end = new_v_encounter.astAb;
SC2.legb.r_in = new_r_encounter.astDa;
SC2.legb.r_end = new_r_encounter.astAb;
SC2.legb.N_rev = sol.Nrev(4);
SC2.legb.TOF = new_tofb*86400/sim.TU;
SC2.legb.Href = sol.Href(:,4);

DT_input = SC2.legb;
tic
[SC2.legb.HS, SC2.legb.RES, SC2.legb.timead, SC2.legb.Xad, SC2.legb.Xpropad, ...
    SC2.legb.Xpropad_DT] = DT_executable(DT_input, sim, data);
el_time.DT_SC2_legb = toc/60; % minutes

SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC2.legb.EPS = marco_elia(SC2.legb,new_r_encounter.astDa,data,sim,SC2.legb.Href,colors);

%% -- if run with less than 100 points
data.n_int = 100;
[SC2.legb] = interpoliamo_porc(SC2.legb,data.n_int);
SC2.legb.v_in = new_v_encounter.astDa;
SC2.legb.v_end = new_v_encounter.astAb;
SC2.legb.r_in = new_r_encounter.astDa;
SC2.legb.r_end = new_r_encounter.astAb;
SC2.legb.N_rev = sol.Nrev(4);
SC2.legb.TOF = new_tofb*86400/sim.TU;
SC2.legb.Href = sol.Href(:,4);
SC2.mjd2000_dep = sol.departure_mjd2000;

SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC2.legb.EPS = marco_elia(SC2.legb,new_r_encounter.astDa,data,sim,SC2.legb.Href,colors);

% -- new mass after DT calculations
SC2.legb.mass_start = SC2.legb.HS.X(1,7);
SC2.legb.mass_end = SC2.legb.HS.X(end,7);
SC2.legb.mass_depleted = SC2.legb.mass_start - SC2.legb.mass_end;
SC2.legb.mass_fraction = SC2.legb.mass_depleted/SC2.legb.mass_start;

% -- coasting on the last asteroid, same coasting time
SC2.mjd2000_dep_asta = SC2.mjd2000_dep+new_tofa+new_cta;
[SC2.coasting.legb] = relative_position_coasting_stuff(SC2.asteroid_b,SC2.mjd2000_dep_asta,...
    new_tofb,new_cta,data.n_int,EphData,sim,colors);

%%
