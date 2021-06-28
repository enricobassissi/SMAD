function [HS, RES, timead, Xad, Xpropad, Xpropad_DT] = DT_executable(DT, sim, data)

% ----- adimensional stuff
DU = sim.DU;
TU = sim.TU;
Mass_start_leg = data.MU;
MU = Mass_start_leg; % it was mdry but we adimensionalise on m_wet
N = data.n_int;

% ---- input related to that leg
v_in = DT.v_in;
v_end = DT.v_end;
r_in = DT.r_in;
r_end = DT.r_end;
N_rev = DT.N_rev;
TOF_adim = DT.TOF;

% ---- run the non linear for this leg
[output] = NL_interpolator_dt( r_in , r_end , v_in , v_end , N_rev , TOF_adim , Mass_start_leg ,sim.PS.Isp , sim);
% ----- interpolate the results within a different discretisation step
output=interpoliamo_insieme_NLI(output,N);

% --- extraction of stuff from NLI
m=output.m;
T=(output.T_magn).*(sim.DU.*1000./sim.TU^2);
r=output.r;
z=output.z;
% s=sqrt(r.^2+z.^2);
vr=output.vr; %enrico deve estrarre il carbone dalle miniere
% vt=output.vt;
vz=output.vz;
% acc_inplane = output.acc_inplane;
% acc_out = output.acc_out;
% acc = output.acc;
TH = output.theta;
gamma = output.gamma;
T_inplane = output.T_inplane.*(sim.DU.*1000./sim.TU^2);
T_outplane = output.T_outplane.*(sim.DU.*1000./sim.TU^2);
theta_dot = output.theta_dot;
time = output.t;

% ------ state vector at nodal points
X = [r, TH, z, vr, theta_dot, vz, m];

% ----- fmincon options
options = optimoptions(@fmincon);
options.Display = 'iter';
options.Algorithm = 'interior-point';
options.MaxIter = 1000;
options.MaxFunEvals = 5e7;   
options.UseParallel = true;
options.ConstraintTolerance = 1e-10;
options.StepTolerance = 1e-13;

XX0 = zeros(1, N*10);
Xad = zeros(N,7);% Xd = Xad;
%adimensionalization of initial conditions
Xad(:,1) = X(:,1);%/DU;
Xad(:,2) = X(:,2); %already adimensional
Xad(:,3) = X(:,3);%/DU;
Xad(:,4) = X(:,4);%/DU*TU; 
Xad(:,5) = X(:,5);%*TU; 
Xad(:,6) = X(:,6);%/DU*TU;
Xad(:,7) = X(:,7)./MU;

% ----- angle of thrust, beta = elevation, alpha = thrust angle in plane
% thr = 1e-3;
% ubeta = atan2((T_outplane < thr).*thr + (T_outplane > thr).*T_outplane, T_inplane);
% ubeta = atan2(T_outplane, T_inplane);
ubeta = asin(T_outplane./ T);
ualpha = gamma;

% --- output parameters results
RES.el = ubeta;
RES.gamma = ualpha;
RES.T_inplane = T_inplane;
RES.T_outplane = T_outplane;
RES.T = T;

% ---- initial and final state
data.xi = Xad(1,:);
data.xf = Xad(end,:);

% ------ propagation with NLI thrust profile
Xpropad = zeros(N, 7);
Xpropad(1,:) = Xad(1,:); %first value equal
timead = time;%/TU;
data.time = timead;

EoMpolarAD=@(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand) EoMpolar(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand,sim);

% RK4 forward integration (in each interval of delta t)
for k = 1:N-1
    xx = Xpropad(k,:);
%         xx = x(k,:);

    hhh = timead(k+1) - timead(k);
    K1 = EoMpolarAD(timead(k), xx, T(k), ualpha(k), ubeta(k), data.muS, data);
    K2 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)), data.muS, data);
    K3 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)),  data.muS, data);
    K4 = EoMpolarAD(timead(k)+ hhh, xx + hhh*K3, T(k+1), ualpha(k+1), ubeta(k+1), data.muS, data);    

%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    Xpropad(k+1,:) = Xpropad(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
end

% ------ definition of the initial guess (sub-optimal NLI solution)
for ii = 1:N
    XX0((ii-1)*10 +1: (ii-1)*10+10) = [Xad(ii,:) T(ii) ualpha(ii) ubeta(ii)];
end

% ---- equality constraints, at initial and final state
Mat_cos=eye(7);
Aeq=[Mat_cos zeros(7,(N*10)-7);zeros(6,(N*10)-10) eye(6) zeros(6,4)];
Beq=[XX0(1:7)';XX0(end-9:end-4)'];  

% ---- definition of upper boundary and lower boundary
[LB, UB] = LBUB(XX0, Xad,  data);

% check to see if any initial state goes over the boundaries from the beginning
% ENRI DOES THINGS: METTI LE BOUNDARIES DINAMICHE CON 1.1 max e min del
% risultato del NLI
% find(XX0 < LB)
% find(XX0 > UB)

% [c, ceq] = EoM(XX0, data, sim);
% nlcon=@(x_hand) EoM(x_hand,data,sim,DT); % thrust directly modulated through power
nlcon=@(x_hand) EoM2(x_hand,data,sim,DT); % thrust modulated throug power and actual thruster datasheet
cost_fun=@(x_hand) DTmethod(x_hand,data);

% ----- call at fmincon
[XSOL,fval,exitflag,fminoutput,lambda,grad,hessian] = fmincon(cost_fun, XX0,[],[],Aeq,Beq,LB,UB,nlcon, options);

% XSOL=fminoutput.bestfeasible.x;
% X = [r, TH, z, vr, theta_dot, vz, m];
XHS = zeros(N, 7);
THS = zeros(1,N);
alphaHS = THS;
betaHS = THS;
% ---- reconstruction of states
% ---- and ri-adimensionality
for ii = 1:length(TH)
    [XHS(ii,:)] = XSOL((ii-1)*10 +1: (ii-1)*10+7);
     XHS(ii,1) = XHS(ii,1) * DU;
     XHS(ii,2) = XHS(ii,2);
     XHS(ii,3) = XHS(ii,3) * DU;
     XHS(ii,4) = XHS(ii,4) * DU/TU;
     XHS(ii,5) = XHS(ii,5) / TU;
     XHS(ii,6) = XHS(ii,6) * DU/TU;
     XHS(ii,7) = XHS(ii,7) * MU;
    [THS(ii)] = XSOL((ii-1)*10+8);
    [alphaHS(ii)] = XSOL((ii-1)*10+9);
    [betaHS(ii)] = XSOL((ii-1)*10+10);
end

rHS = XHS(:,1);

% ---- output parameters Happy Solution
HS.X = XHS;
HS.r = rHS;
HS.beta = betaHS'; % because they would be row vectors
HS.alpha = alphaHS';
HS.T = THS';

% ---- check physical validity, propagation with DT results
Xpropad_DT(1,:) = Xad(1,:);
% RK4 forward integration (in each interval of delta t)
for k = 1:N-1
    xx = Xpropad_DT(k,:);
%         xx = x(k,:);

    hhh = timead(k+1) - timead(k);
    K1 = EoMpolarAD(timead(k), xx, THS(k), alphaHS(k), betaHS(k), data.muS, data);
    K2 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K1, 0.5*(THS(k+1) + THS(k)), 0.5*(alphaHS(k+1) + alphaHS(k)), 0.5*(betaHS(k)+betaHS(k+1)), data.muS, data);
    K3 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K2, 0.5*(THS(k+1) + THS(k)), 0.5*(alphaHS(k+1) + alphaHS(k)), 0.5*(betaHS(k)+betaHS(k+1)),  data.muS, data);
    K4 = EoMpolarAD(timead(k)+ hhh, xx + hhh*K3, THS(k+1), alphaHS(k+1), betaHS(k+1), data.muS, data);    

%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    Xpropad_DT(k+1,:) = Xpropad_DT(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
end

end