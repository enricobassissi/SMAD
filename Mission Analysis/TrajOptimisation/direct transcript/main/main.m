%% work environment setup
clear all
str_path=split(pwd, 'TrajOptimisation\direct transcript\main');
util_path=string(str_path(1))+'Utils';
addpath(genpath(util_path));
str_path_1=split(pwd, 'main');
imp_path=string(str_path_1(1))+'functions';
addpath(genpath(imp_path));
%% Default options
set(0, 'DefaultTextFontSize', 12) % modify it if too small
set(0, 'DefaultAxesFontSize', 11) % modify it if too small
set(0, 'DefaultLegendFontSize', 12) % modify it if too small
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 1.2)

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
      
%% LOAD INIT GUESS DATA
load('ws_2RL_all_indietro_moo2.mat')
sim.M1=output.m_SC1(1);
clear data output

%% DIRECT TRANSCRIPTION
N=50;
v_launcher = sol.v_inf_magn/sim.DU*sim.TU*[cos(sol.el)*cos(sol.az); cos(sol.el)*sin(sol.az); sin(sol.el)];
v_dep = v_encounter.EA + v_launcher;  %if parabolic escape (v_extra = 0)
[output] = NL_interpolator_dt( r_encounter.EA , r_encounter.astA1 , v_dep , v_encounter.astA1 , sol.Nrev(1) , sol.TOF1_ADIM ,sim.M1 ,sim.PS.Isp ,sim);
output=interpoliamo_insieme_NLI(output,N);
data.n_int = N; % discretisation points
data.Tmax = 0.025; % max thrust available
href = output.Href; % vector normal to the transfer orbit plane
TOF = sol.TOF1; % time of flight, dimensional [days]
t0 = sol.departure_mjd2000; %to be checked, maybe 0
r1vers = r_encounter.EA./(norm(r_encounter.EA));

m=output.m;
T=(output.T_magn).*(sim.DU.*1000./sim.TU^2);
r=output.r;
z=output.z;
s=sqrt(r.^2+z.^2);
vr=output.vr; %enrico deve estrarre il carbone dalle miniere
vt=output.vt;
vz=output.vz;
acc_inplane = output.acc_inplane;
acc_out = output.acc_out;
acc = output.acc;
TH = output.theta;
gamma = output.gamma;
T_inplane = output.T_inplane.*(sim.DU.*1000./sim.TU^2);
T_outplane = output.T_outplane.*(sim.DU.*1000./sim.TU^2);
theta_dot = output.theta_dot;
time = output.t;
% for i=1:length(T)
% t_onoff=0
% end;
% ENRI: L, gamma1, gamma2, v1perp, v2perp, v1tra, v2tra, vnorm, dmdt, TOFr || non li usa mai!

% [ m, T, r, z, s, vr, vt, vz, acc_inplane, acc_out, acc, TH, L, gamma1, gamma2, gamma, v1perp, v2perp, v1tra, v2tra, vnorm, dmdt, T_inplane, T_outplane, theta_dot, time, TOFr] = ...
%     Conway(TOF3, N_rev3, q3, r1norm_3, r2norm_3, r1vers_3, r2vers_3, hvers_3, hh_3, v1_3, v2_3, muS, data);
X = [r, TH, z, vr, theta_dot, vz, m];

% X = [r', TH', z', vr', theta_dot', vz', m'];
options = optimset( 'Algorithm','interior-point', ...
                        'Display','iter', ... %          
                        'MaxIter',250, ...
                        'LargeScale','on', ...
                        'MaxFunEvals', 500000);
%                         
options.ConstraintTolerance=1e-10;
DU = sim.DU;
TU = sim.TU;
MU = sim.M1; % ENRI: it was mdry but we adimensionalise on mwet
% 
% DU = 1;
% TU = 1;
% MU = 1;

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

% thr = 1e-3;
% ubeta = atan2((T_outplane < thr).*thr + (T_outplane > thr).*T_outplane, T_inplane);
% ubeta = atan2(T_outplane, T_inplane);
ubeta = asin(T_outplane./ T);
ualpha = gamma;

data.xi = Xad(1,:);
data.xf = Xad(end,:);
data.muS = sim.mu_dim;   
data.MU = sim.M1; % ENRI: same thing on mass adimensionalisation
data.Isp = sim.PS.Isp; % ENRI: our isp

%%%ode better? #TODO
% acceleration=@acc_stepwise(
% opts = odeset('Reltol',1e-12,'AbsTol',1e-14,'Stats','on');
% [T,state_sol]=ode113(@(t,X) EoMpolar(~, x, T,ualpha(t), beta(t), ~, data,sim), output.t, x_,opts);
% %%%

Xpropad = zeros(N, 7);
Xpropad(1,:) = Xad(1,:); %first value equal
timead = time;%/TU;
data.time = timead;

EoMpolarAD=@(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand) EoMpolar(t_hand, x_hand, T_hand, alpha_hand, beta_hand, muS_hand, data_hand,sim);

    %RK4 forward integration (in each interval of delta t)
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
% Xd(:,1) = Xpropad(:,1);%*DU;
% Xd(:,2) = Xpropad(:,2);%; %already adimensional
% Xd(:,3) = Xpropad(:,3);%*DU;
% Xd(:,4) = Xpropad(:,4);%*DU/TU; 
% Xd(:,5) = Xpropad(:,5);%/TU; 
% Xd(:,6) = Xpropad(:,6);%*DU/TU;
% Xd(:,7) = Xpropad(:,7)/MU;

figure()
subplot(7,1,1), plot(Xad(:,1)), hold on, plot(Xpropad(:,1))
subplot(7,1,2), plot(Xad(:,2)), hold on, plot(Xpropad(:,2))
subplot(7,1,3), plot(Xad(:,3)), hold on, plot(Xpropad(:,3))
subplot(7,1,4), plot(Xad(:,4)), hold on, plot(Xpropad(:,4))
subplot(7,1,5), plot(Xad(:,5)), hold on, plot(Xpropad(:,5))
subplot(7,1,6), plot(Xad(:,6)), hold on, plot(Xpropad(:,6))
subplot(7,1,7), plot(Xad(:,7)), hold on, plot(Xpropad(:,7)), legend('ad','prop')

%already adimensional
data.xi = Xad(1,:);
data.xf = Xad(end,:);

% Adimensionalization of variables
figure(),
subplot(5,1,1), plot(rad2deg(ubeta)), title('Beta'), grid on
subplot(5,1,2), plot(rad2deg(gamma)), title('Gamma'), grid on
subplot(5,1,3), plot(T_inplane), title('T inplane'), grid on
subplot(5,1,4), plot(T_outplane), title('T outplane'), grid on
subplot(5,1,5), plot(T), title('T'), grid on
 
%definition of the initial guess (sub-optimal conway solution)
for ii = 1:N
    XX0((ii-1)*10 +1: (ii-1)*10+10) = [Xad(ii,:) T(ii) ualpha(ii) ubeta(ii)];
end
%boundaries NOT IN USE
Mat_cos=eye(7);%-[0 0 0 0 0 ; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 0;0 0 0 0 0; ];
Aeq=[Mat_cos zeros(7,(N*10)-7);zeros(6,(N*10)-10) eye(6) zeros(6,4)];
Beq=[XX0(1:7)';XX0(end-9:end-4)'];       
%%
% kepMtry = uplanet(t0 + TOF,4);
% rMtry = kep2car2(kepMtry, muS);
% R = refplane2car( X(end,1), X(end,3),  X(end,1)*X(end,5), X(end,4), X(end,6), X(end,2), r1vers, href);
% rMtry - R
%definition of upper boundary and lower boundary
[LB, UB] = LBUB(XX0, Xad,  data);
find(XX0 < LB)
find(XX0 > UB)
[c, ceq] = EoM(XX0, data, sim);
nlcon=@(x_hand) EoM(x_hand,data,sim);
cost_fun=@(x_hand) DTmethod(x_hand,data);
[XSOL,fval,exitflag,fminoutput,lambda,grad,hessian] = fmincon(cost_fun, XX0,[],[],Aeq,Beq,LB,UB,nlcon, options);
%%

%XSOL=fminoutput.bestfeasible.x;
XHS = zeros(N, 7);
THS = zeros(1,N);
alphaHS = THS;
betaHS = THS;
%reconstruction of states
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
figure()
sgtitle('Control allocation')
subplot(3,1,1), plot(timead, THS), hold on, plot(timead, T), hold on, title('T')
subplot(3,1,2), plot(timead, rad2deg(alphaHS)), hold on, plot(timead, rad2deg(ualpha)), hold on, title('alpha')
subplot(3,1,3), plot(timead, rad2deg(betaHS)), hold on, plot(timead, rad2deg(ubeta)), hold on, title('beta')
% figure()
%plot(timead, rHS), hold on, plot(timead, r), hold on,
%% CHECK PHYSICS VALIDITY
clear Xpropad
Xpropad(1,:) = Xad(1,:);
    %RK4 forward integration (in each interval of delta t)
    for k = 1:N-1
        xx = Xpropad(k,:);
%         xx = x(k,:);

        hhh = timead(k+1) - timead(k);
        K1 = EoMpolarAD(timead(k), xx, THS(k), alphaHS(k), betaHS(k), data.muS, data);
        K2 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K1, 0.5*(THS(k+1) + THS(k)), 0.5*(alphaHS(k+1) + alphaHS(k)), 0.5*(betaHS(k)+betaHS(k+1)), data.muS, data);
        K3 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K2, 0.5*(THS(k+1) + THS(k)), 0.5*(alphaHS(k+1) + alphaHS(k)), 0.5*(betaHS(k)+betaHS(k+1)),  data.muS, data);
        K4 = EoMpolarAD(timead(k)+ hhh, xx + hhh*K3, THS(k+1), alphaHS(k+1), betaHS(k+1), data.muS, data);    

%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
        Xpropad(k+1,:) = Xpropad(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
% Xd(:,1) = Xpropad(:,1);%*DU;
% Xd(:,2) = Xpropad(:,2);%; %already adimensional
% Xd(:,3) = Xpropad(:,3);%*DU;
% Xd(:,4) = Xpropad(:,4);%*DU/TU; 
% Xd(:,5) = Xpropad(:,5);%/TU; 
% Xd(:,6) = Xpropad(:,6);%*DU/TU;
% Xd(:,7) = Xpropad(:,7)/MU;

figure()
subplot(7,1,1), plot(Xad(:,1)), hold on, plot(Xpropad(:,1))
subplot(7,1,2), plot(Xad(:,2)), hold on, plot(Xpropad(:,2))
subplot(7,1,3), plot(Xad(:,3)), hold on, plot(Xpropad(:,3))
subplot(7,1,4), plot(Xad(:,4)), hold on, plot(Xpropad(:,4))
subplot(7,1,5), plot(Xad(:,5)), hold on, plot(Xpropad(:,5))
subplot(7,1,6), plot(Xad(:,6)), hold on, plot(Xpropad(:,6))
subplot(7,1,7), hold on, plot(Xpropad(:,7)), legend('conway','propagate dt')

%%
function J = DTmethod(X, data)

%     N =data.n_int;
    %T = m;
    %time = data.time;
% alpha = m; beta = m;
    m(1)=X(7);
    m(2)=X(end-3);
%     
%     J = 0;
%     for k = 2:N-1
%         %integration of cost function
%         dt = time(k+1) - time(k);
%         J = J +  dt/6 * (T(k-1) + 4*T(k) + T(k + 1));
%     end
%     J = J + T(end)*dt;
%     T_inplane = a_in .*m * 1000;
%     T_outplane = a_out .*m * 1000;
%     T = (T_inplane.^2 + T_outplane.^2).^0.5;
    J = (m(1) - m(2))*data.MU;

end


function [lb, ub] = LBUB(XX0, X,  data)
%X0 vector of propagated ODE45
    N = data.n_int;
    T_lb = 0;
    T_ub = data.Tmax;
    alpha_lb = deg2rad(-360);%-pi;
    alpha_ub = deg2rad(360);%pi;
    beta_lb = deg2rad(-360);%-pi/4;
    beta_ub = deg2rad(360);%pi/4;
    
    %limits definition
    ub = zeros(1,10*N);
    lb = zeros(1,10*N);
        
    %r
    r_ub = 1.5;
    r_lb = 0.8;

    %theta
    th_ub = 20;
    th_lb = 0;
    %z
    z_ub =  0.05;
    z_lb =  -0.05;
    %vr
    vr_ub = 0.15;
    vr_lb = -0.2;
    %theta_dot
    thd_ub = Inf;
    thd_lb = -Inf;
    %vz
    vz_ub = 0.03;
    vz_lb = -0.03;
    
    %
    m_ub = X(1,7)*1.1;
    m_lb = 0.5;

    %INITIAL CONDITION
    ub(1:10) = XX0(1:10); 
    ub(8) = T_ub; ub(9) = alpha_ub; ub(10) = beta_ub;
    ub(7) = m_ub;
    lb(1:10) = XX0(1:10); 
    lb(7) = m_lb;
    lb(8) = 0; lb(9) = alpha_lb; lb(10) = beta_lb;
    
    
    %PATH CONSTRAINTS
    ub(11:10:((N-1)*10)) = r_ub;
    ub(12:10:((N-1)*10)) = th_ub;
    ub(13:10:((N-1)*10)) = z_ub;
    ub(14:10:((N-1)*10)) = vr_ub;
    ub(15:10:((N-1)*10)) = thd_ub;
    ub(16:10:((N-1)*10)) = vz_ub;
    ub(17:10:((N-1)*10)) = m_ub;
    ub(18:10:((N-1)*10)) = T_ub;
    ub(19:10:((N-1)*10)) = alpha_ub;
    ub(20:10:((N-1)*10)) = beta_ub;

    lb(11:10:((N-1)*10)) = r_lb;
    lb(12:10:((N-1)*10)) = th_lb;
    lb(13:10:((N-1)*10)) = z_lb;
    lb(14:10:((N-1)*10)) = vr_lb;
    lb(15:10:((N-1)*10)) = thd_lb;
    lb(16:10:((N-1)*10)) = vz_lb;
    lb(17:10:((N-1)*10)) = m_lb;%data.Mdry;
    lb(18:10:((N-1)*10)) = T_lb;
    lb(19:10:((N-1)*10)) = alpha_lb;
    lb(20:10:((N-1)*10)) = beta_lb;

    %FINAL CONDITION
    ub(end-9:end) = XX0(end-9:end);
    ub(end-2)=T_ub; ub(end-1)=alpha_ub; ub(end)=beta_ub;
    ub(end-3) = m_ub;
    lb(end-9:end) = XX0(end-9:end);
    lb(end-2)=0; lb(end-1)=alpha_lb; lb(end)=beta_lb;
    lb(end-3) = m_lb;
end

%MY ACC STEPWISE MM
function out = acc_stepwise(u,t,t_vec)
k=0;
for i=1:(length(t_vec)-1)
if and(t_vec(i)<=t,t<t_vec(i+1))
    k=i;
end  
end
if k==0
    out=0;
else
    out=(u(k));
end
end