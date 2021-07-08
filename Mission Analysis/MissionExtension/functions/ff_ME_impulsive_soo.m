function obj_fun = ff_ME_impulsive_soo(x, data, sim)
% setting the input times
MJD01 = x(1)/86400*sim.TU;
TOF = x(2)/86400*sim.TU;
MJDF1 = MJD01 + TOF;

% Computing position and velocity of the asteroids in those moments
% departure from 1st ast
state_ast_1 = interp1(data.time_eval,data.horizons_ast1,MJD01); % AU, AU/day
r1 = state_ast_1(1:3); % au
v1 = state_ast_1(4:6)/86400*sim.TU; % adim

% arrival at 2nd asteroid 
state_ast_2 = interp1(data.time_eval,data.horizons_ast2,MJDF1); % AU, AU/day
r2 = state_ast_2(1:3);
v2 = state_ast_2(4:6)/86400*sim.TU; % adim

[dv,~,~]=lambert_solver_rendezvous(r1,r2,v1,v2,TOF*86400/sim.TU,sim.mu); % everything adim... AU, adim vel, adim mu
% dvtot = dv*sim.DU/sim.TU; % ri adimensionalise -> km/s
obj_fun = dv;

end

