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
%%
load('soo_ps_flyby_3-33_1dayeachpointTraj.mat');%load('soo_ps_flyby_3-33.mat');
sub_sol=zeros(4,2);
fval=zeros(4,1);

%% PS param
options = optimoptions('particleswarm');
options.HybridFcn = @fmincon;
options.SwarmSize = 1000; % Default is min(100,10*nvars),
options.MaxIterations = 100; %  Default is 200*nvars
options.MaxStallIterations = 50; % Default 20
options.Display = 'iter';
options.FunctionTolerance = 1e-6;

%% landing asteroid 1
asteroid=sol.ast_1;
depdate=sol.MJD0;
obj_fun=@(state) ff_balistic_landing(state,sol,asteroid,muSun,depdate,data);
lb=[1;x(2)];
ub=[x(2)-1;x(2)+20];
% opts=optimoptions('ga','MaxGeneration',ngen,'FunctionTolerance',10e-10);
% [sub_sol(1,:),fval(1)]=ga(obj_fun,2,[],[],[],[],lb,ub,[],opts)


[sub_sol(1,:),fval(1),exitFlag,Output] = particleswarm(obj_fun,2,lb,ub,options);
%% landing asteroid 2
asteroid=sol.ast_2;
depdate=sol.MJD0+x(2);
obj_fun=@(state) ff_balistic_landing(state,sol,asteroid,muSun,depdate,data);
lb=[x(3)-20;x(3)];
ub=[x(3)-1;x(3)+20];
[sub_sol(2,:),fval(2),exitFlag,Output] = particleswarm(obj_fun,2,lb,ub,options);
%% landing asteroid 3
asteroid=sol.ast_3;
depdate=sol.MJD0+x(2)+x(3);
obj_fun=@(state) ff_balistic_landing(state,sol,asteroid,muSun,depdate,data);
lb=[x(4)-20;x(4)];
ub=[x(4)-1;x(4)+20];
[sub_sol(3,:),fval(3),exitFlag,Output] = particleswarm(obj_fun,2,lb,ub,options);

%% landing asteroid 4
asteroid=sol.ast_4;
depdate=sol.MJD0+x(2)+x(3)+x(4);
obj_fun=@(state) ff_balistic_landing(state,sol,asteroid,muSun,depdate,data);
lb=[x(5)-20;x(5)];
ub=[x(5)-1;x(5)+20];
[sub_sol(4,:),fval(4),exitFlag,Output] = particleswarm(obj_fun,2,lb,ub,options);
