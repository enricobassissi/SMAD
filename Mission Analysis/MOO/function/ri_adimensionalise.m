function x = ri_adimensionalise(x,sim)

x(:,1) = x(:,1)*sim.TU/(3600*24); % MJD01
x(:,2) = x(:,2)*sim.TU/(3600*24); % TOF1
x(:,3) = x(:,3)*sim.DU/sim.TU; % v_inf_magn
% alpha and beta were already fine
x(:,6) = x(:,6)*sim.TU/(3600*24); % Buffer time 1
x(:,7) = x(:,7)*sim.TU/(3600*24); % TOF2
sim.bound.TOF2_min = 50*3600*24/sim.TU; % days
sim.bound.TOF2_max = 3*365*3600*24/sim.TU; % days
% Matrix of permutations is already fine
x(:,9) = x(:,9)*sim.TU/(3600*24); % Buffer time 2
x(:,10) = x(:,10)*sim.TU/(3600*24);  TOF3
x(:,11) = x(:,11)*sim.TU/(3600*24); % Buffer time 3 
x(:,12) = x(:,12)*sim.TU/(3600*24); % TOF4

end