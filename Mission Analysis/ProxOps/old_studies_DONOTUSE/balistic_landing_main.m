clear all
%% add path of functions and python stuff
str_path=split(pwd, 'ProxOps\asteroid relative dynamics');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
py_path=string(str_path(1))+'PyInterface\NEO_API_py';
addpath(genpath(py_path));
neoeph_path=string(str_path(1))+'NeoEph';
addpath(genpath(neoeph_path));
trajopt_path=string(str_path(1))+'TrajOptimisation';
addpath(genpath(trajopt_path));
str_path=split(pwd, 'main');
imp_path=string(str_path(1));
addpath(genpath(imp_path));

%% Call to NASA JPL Horizons to get Asteroid's Ephemerides
% Import module of Python
try 
    module = py.importlib.import_module('neo_api_function');
catch
    copyfile(py_path+'\neo_api_function.py', pwd, 'f'); 
    module = py.importlib.import_module('neo_api_function');
end

%% Asteroids
AU = astroConstants(2);
muSun = astroConstants(4);
asteroid_1 = "2012SY49";%"2020UE";
asteroid_2 = "2005WG57";%"2012QD8";
asteroid_3 = "2012QD8";%"2012SY49";
asteroid_4 = "2008XU2";
%["2012SY49","2005WG57","2012QD8","2008XU2"]
%["2012SY49","2006HX57","2012BY1","2020UE"]
% data extraction section
data.asteroid_names = [asteroid_1,asteroid_2,asteroid_3,asteroid_4];
[data.y_interp_ft, data.t_vector] = find_eph_neo(data.asteroid_names);

x=[9.081585446108058e+03,8.691662798408688e+02,5.181098966713850e+02,5.309387183646369e+02,6.911875175294085e+02,1.280338755302122e+03]
%x=[9.6127e+03, 359.5569, 440.6470, 638.0827, 727.2520]
%[9.492266228410985e+03,2.675577692802069e+02,1.125558390299206e+03,1.115318519817115e+03,4.888829584204003e+02,1048,0,0,0,0]


Sun=body;
Sun=Sun.Sun;

SC=body;
SC=SC.Spacecraft;

%%
if length(x)==5
% Departure from Earth
[kep_EA,ksun] = uplanet(x(1), 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
% passage at 1st ast
[kep_ast_1] = uNEO2(x(1)+x(2),asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
[~,~,~,~,V0_sc,~,~,~] = lambertMR(rEA,r1,x(2)*3600*24,ksun,0,0,0,0);
obj_fun = @(state_vector) eval_balistic_landing(state_vector,data,x(1),asteroid_1,rEA,V0_sc,Sun,SC);
%% options=optimoptions('ga','MaxGenerations',400,'MaxTime',3600) 
[sol,fval]=ga(obj_fun,2,[],[],[],[],[30 x(2)+1],[x(2)-1 x(2)+4])
%%
mothership_orbit=SC.propagate_around(Sun,rEA,V0_sc,sol(2)*3600*24);
plot3(mothership_orbit.r(:,1),mothership_orbit.r(:,2),mothership_orbit.r(:,3),'-')
hold on
plot3(r1(1),r1(2),r1(3),'x')
for i=1:length(mothership_orbit.T)
    if and(mothership_orbit.T(i)<sol(1)*3600*24,mothership_orbit.T(i+1)>=sol(1)*3600*24)
        break
    end
end
[~,~,~,~,VI,VF,~,~] = lambertMR(mothership_orbit.r(i,:),r1,sol(2)*3600*24,Sun.mu,0,0,0,0);
sc_orbit=SC.propagate_around(Sun,mothership_orbit.r(i,:),VI,sol(2)*3600*24);
plot3(sc_orbit.r(:,1),sc_orbit.r(:,2),sc_orbit.r(:,3))
hold on
else
    % Departure from Earth
[kep_EA,ksun] = uplanet(x(1)+x(2), 3);
[rEA, vEA] = sv_from_coe(kep_EA,ksun); % km, km/s
% passage at 1st ast
[kep_ast_1] = uNEO2(x(1)+x(2)+x(3),asteroid_1,data); % [km,-,rad,rad,rad,wrapped rad]
[r1, v1] = sv_from_coe(kep_ast_1,ksun); % km, km/s
[~,~,~,~,V0_sc,~,~,~] = lambertMR(rEA,r1,x(3)*3600*24,ksun,0,0,0,0);
obj_fun = @(state_vector) eval_balistic_landing(state_vector,data,x(1)+x(2),asteroid_1,rEA,V0_sc,Sun,SC);
%options=optimoptions('ga','MaxGenerations',400,'MaxTime',3600) 
opts=optimoptions('ga', 'PopulationSize',40,'MaxGenerations',200);
[sol,fval]=ga(obj_fun,2,[],[],[],[],[1; x(3)],[x(3)-1; x(3)+30],[],opts);
end
%%
mothership_orbit=SC.propagate_around(Sun,rEA,V0_sc,sol(2)*3600*24);
plot3(mothership_orbit.r(:,1),mothership_orbit.r(:,2),mothership_orbit.r(:,3),'-')
hold on
plot3(r1(1),r1(2),r1(3),'x')
for i=1:length(mothership_orbit.T)
    if and(mothership_orbit.T(i)<sol(1)*3600*24,mothership_orbit.T(i+1)>=sol(1)*3600*24)
        break
    end
end
[~,~,~,~,VI,VF,~,~] = lambertMR(mothership_orbit.r(i,:),r1,(sol(2)-sol(1))*3600*24,Sun.mu,0,0,0,0);
sc_orbit=SC.propagate_around(Sun,mothership_orbit.r(i,:),VI,sol(2)*3600*24);
plot3(sc_orbit.r(:,1),sc_orbit.r(:,2),sc_orbit.r(:,3))
hold on