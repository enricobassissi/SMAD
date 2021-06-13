%% work environment setup
clear all
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
%% load workspace
%load Conway results workspace
load('topputo_fittizio.mat')
%load precomputed Gradient for Nonlinear constraint on residual 
% load('d_stored.mat')
load('d_stored.mat')
%transform precomputed Gradient from Simbolic Toolbox obj to handle function
d_residual=matlabFunction(d_residual);
% clear ambiguos variables
clear lb ub Aeq Beq A B nlcon fval Fval n N r u v options sim
%% simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7           ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 500;                    % number of computational nodes
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.PS.Isp = 3000/sim.TU;
N=sim.n_sol;
%% extract orbit data from Conway
%init vectors
r=zeros(sim.n_sol,2);
v=zeros(sim.n_sol,2);
%extract r as [teta r]
for i=1:sim.n_sol
r(i,:)=[output.l(i) output.r(i)];
end
%extract velocities
for i=1:sim.n_sol
    %vr
    v(i,2)=output.v_r(i);  
    %vt
    v(i,1)=output.w(i)*output.r(i);
    %ut
    u(i)=((sim.TU^2)/(1000*sim.DU))*output.u(i,1)/output.m(i);%((sum(u_cart(i,1:2).*R_norm_car)/(norm(R_norm_car)^2))*R_norm_car);  
end
u(1)=0;
u(sim.n_sol)=0;


%inserted time discretization
tin=0; %adim
tfin=output.t(N); %adim
time_vec=output.t;%adim

%maximum adim tang acceleration
ac=max(u)*1.01;

%assembly Conway orbit data to state vector [teta1 vt1 u1 r1 vr1 teta2
%vt2...]
init_guess=[];
for i=1:N
    init_guess=[init_guess; r(i,1);v(i,1);u(i);r(i,2);v(i,2)];
end

%boundaries conditions
initial_states=[r(1,1),v(1,1),0,r(1,2),v(1,2)]';
final_states=[r(N,1),v(N,1),0,r(N,2),v(N,2)]';
Mat_cos=eye(5);
Aeq=[Mat_cos zeros(5,(N*5)-5);zeros(5,(N*5)-5) Mat_cos];
Beq=[initial_states;final_states];
%lower and upper bounds
max_V=2;
for i=1:N
    n=(i-1)*5;
    
    lb(n+1)=init_guess(n+1)-2*pi;
    lb(n+4)=init_guess(n+4)-0.5; 
    
    lb(n+2)=0; 
    lb(n+5)=-max_V; 
    
    lb(n+3)=0; 
   
    ub(n+1)=init_guess(n+1)+2*pi;
    ub(n+4)=init_guess(n+4)+0.5;
    
    ub(n+2)=max_V;
    ub(n+5)=max_V;
    
    ub(n+3)=ac; 
    
end


%cost fumction (mass ratio)
cost_fun=@(x) thrust_integrator(x,N,time_vec,output.m(1),tin,tfin,sim);

%non linear constraint (residuals of the Hermite Simpson collocation)
nlcon=@(x) non_lin_con(x,N,time_vec,d_residual);
% nlcon=@(x) non_lin_Lobatto(x,N,time_vec);

A=[]; b=[];

%options
options = optimoptions('fmincon');
options.MaxFunctionEvaluations =1e+06;
options.MaxIter=500;
options.Display='iter';
options.UseParallel=true;

% options.CheckGradients=true;
options.SpecifyConstraintGradient=true;
options.ConstraintTolerance=1e-6;
options.StepTolerance=1e-12;


%% solve
[x_sol,fval,eflag,fminoutput]=fmincon(cost_fun,init_guess,A,b,Aeq,Beq,lb,ub,nlcon,options);
%% plot
for i=1:N-1
    n=(i-1)*5;
    x_k(i)=x_sol(n+1); 
    y_k(i)=x_sol(n+4); 
    [x_k(i),y_k(i)]=pol2cart(x_k(i),y_k(i));
    vx_k(i)=x_sol(n+2); 
    vy_k(i)=x_sol(n+5); 
    ux_k(i)=x_sol(n+3);
end
figure(1)
plot(x_k,y_k,'Color',colors(1,:));
figure(2)
bar(ux_k);