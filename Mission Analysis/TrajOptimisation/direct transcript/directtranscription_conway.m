
clear all
% load('R3dim.mat')
% load('V3dim.mat')
% load('Tglobal.mat')
% %insert here the r history vector spanned uniformly
% r=R3_dim;%[0 0 1; 0 0 2];
% %insert here the v history vector spanned uniformly
% v=V3_dim;%[0 0 1; 0 0 0];
% load('R3400.mat')
% load('V3400.mat')
% load('Tglobal400.mat')
% load('outputm400.mat')
% load('lo_sballo.mat');
% output.Thrust=output.u;
% output.theta=output.l
% =output.l
% load('mars_lt_fordirect');
% load('ws_1_ast_0-09.mat')
load('topputo_fittizio.mat')
clear lb ub Aeq Beq A B nlcon fval Fval n N r u v options sim
%% simulation parameters
sim.mu_dim    = 132712440018          ; % actractor parameter [km^3 s^-2]
sim.DU        = 149597870.7           ; % distance unit [km]
sim.TU        = (sim.DU^3/sim.mu_dim )^0.5; % time unit [s]
sim.mu        = 1;                      % non-dimensional attractor parameter [DU^3/TU^2]
sim.n_sol     = 500;                    % number of computational nodes
sim.g0 = 9.81*(sim.TU^2/(1000*sim.DU)); % non-dimensional g0
sim.PS.Isp = 3000/sim.TU;
%%
%insert here the r history vector spanned uniformly
r=zeros(sim.n_sol,2);
v=zeros(sim.n_sol,2);
%u_cart=((sim.TU^2)/(1000*sim.DU))*Tlocal(:,1:2)./output.m;%[0 0 -1; 0 0 0]; adim acc
% u_=[R3(1,1),R3(1,2),R3(1,3)]';
for i=1:sim.n_sol %1 az 2 el 3 ro
%     [r(i,1),r(i,2)]=cart2pol(R3(i,1),R3(i,2));
% v_(:)=[R3(i,1),R3(i,2),R3(i,3)]';
% theta_=atan2(norm(cross(u_,v_)),dot(u_,v_));%vecangle(u_,v_);
%     if vecangle(u_,v_)<0
%         r(i,:)=[pi+(pi-theta_);norm(v_)];  % cart2pol(R3(i,1),R3(i,2));
%     else
%         r(i,:)=[theta_;norm(v_)];
%     end
% end
r(i,:)=[output.l(i);output.r(i)];
end
r(:,1)=unwrap(r(:,1));

%UNCOMMENT THIS
% for i=1:sim.n_sol
%     V_norm(i)=output.vin(i);
%     %R_car=[r3(i,1),r3(i,2)];
%     %R_norm_car=[-r3(i,2),r3(i,1)];
%     %vt
%     v(i,1)=V_norm(i)*cos(output.Thrust(i,2));%norm((sum(V_car.*R_norm_car)/(norm(R_norm_car)^2))*R_norm_car);  
%     %Vr
%     v(i,2)=V_norm(i)*sin(output.Thrust(i,2));%norm((sum(V_car.*R_car)/(norm(R_car)^2))*R_car);
%     %vz
%     %v(i,3)=direct_data.V3(i,3);
%     %ut
%     u(i)=((sim.TU^2)/(1000*sim.DU))*output.Thrust(i,1)/output.m(i);%((sum(u_cart(i,1:2).*R_norm_car)/(norm(R_norm_car)^2))*R_norm_car);  
%     %ur
%     %u_not_considered(i)=norm((sum(u_cart(i,i:2).*R_car)/(norm(R_car)^2))*R_car);
%     %uz
%     %u(i,3)=u_cart(i,3);
% end
%%END UNCOMMENT


for i=1:sim.n_sol
%     V_car(:)=V3(i,1:2);
%     R_car=[r3(i,1),r3(i,2)];
%     R_norm_car=[-r3(i,2),r3(i,1)];
    
%     V_car(:)=V3(i,1:3);
%     R_car=R3(i,:);
%     R_norm_car=[-r3(i,2),r3(i,1)];
    %Vr
    v(i,2)=output.v_r(i);  %norm((sum(V_car.*R_car)/(norm(R_car)^2))*R_car);%  
    %vt
    v(i,1)=output.w(i)*output.r(i);%+sqrt(1/r(1,2));%sqrt(norm(V_car)^2-v(i,2)^2); norm((sum(V_car.*R_norm_car)/(norm(R_norm_car)^2))*R_norm_car);  %V_norm(i)*cos(output.Thrust(i,2));

    %vz
    %v(i,3)=direct_data.V3(i,3);
    %ut
    u(i)=((sim.TU^2)/(1000*sim.DU))*output.u(i,1)/output.m(i);%((sum(u_cart(i,1:2).*R_norm_car)/(norm(R_norm_car)^2))*R_norm_car);  
    %ur
    %u_not_considered(i)=norm((sum(u_cart(i,1:2).*R_car)/(norm(R_car)^2))*R_car);
    %uz
    %u(i,3)=u_cart(i,3);
end



%insert here the control history vector spanned uniformly


%insert initial time, final time and time resolution for the previous
%inserted data
N=sim.n_sol;
tin=0;
tfin=output.t(N);
time_vec=output.t;
%insert max thrust
%ac=0.24./(100*1000); %GUARDA TODO SOTTO
ac=((sim.TU^2)/(1000*sim.DU))*max(u)/1000;%NEWTON to adim acc for a default 1000kg mass
%compose the initial guess vector
init_guess=[];
for i=1:N
    %init_guess=[init_guess; r(i,:)'; v(i,:)'; u(i,:)'];
    init_guess=[init_guess; r(i,1);v(i,1);u(i);r(i,2);v(i,2)];
end
%time interval duration
h=(tfin-tin)/(N-1);
%boundaries conditions
initial_states=[r(1,1),v(1,1),0,r(1,2),v(1,2)]';
final_states=[r(N,1),v(N,1),0,r(N,2),v(N,2)]';
%assembling the problem for the NLP solver
cost_fun=@(x) thrust_integrator(x,N,time_vec,output.m(1),tin,tfin,sim);
nlcon=@(x) non_lin_con(x,N,time_vec);
% nlcon=@(x) non_lin_Lobatto(x,N,time_vec);
Mat_cos=eye(5);%-[0 0 0 0 0 ; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 0;0 0 0 0 0; ];
Aeq=[Mat_cos zeros(5,(N*5)-5);zeros(5,(N*5)-5) Mat_cos];
Beq=[initial_states;final_states];
A=[]; b=[]; lb=[]; ub=[];
options = optimoptions('fmincon');
options.MaxFunctionEvaluations =5e+04;
options.Display='iter';
options.UseParallel=true;
options.MaxIter=100;
options.ConstraintTolerance=1;
% options.OptimalityTolerance=1e-12;
% options.FunctionTolerance=1e-12;
options.StepTolerance=1e-13;
max_V=2*2;
for i=1:N
%     n=(i-1)*9;
%     lb(n+1)=-pi; 
%     lb(n+4)=-pi/2; 
%     lb(n+7)=0; 
%     lb(n+2)=-max_V; 
%     lb(n+5)=-max_V; 
%     lb(n+8)=-max_V;
%     lb(n+3)=-ac; 
%     lb(n+6)=-ac; 
%     lb(n+9)=-ac; 
%     ub(n+1)=pi; 
%     ub(n+4)=pi/2; 
%     ub(n+7)=2; 
%     ub(n+2)=max_V; 
%     ub(n+5)=max_V; 
%     ub(n+8)=max_V;
%     ub(n+3)=ac; 
%     ub(n+6)=ac; 
%     ub(n+9)=ac;
    n=(i-1)*5;
    lb(n+1)=0;%init_guess(n+1)-4*pi;%0.26/30; 
    lb(n+4)=0; 
    
    lb(n+2)=0; 
    lb(n+5)=-max_V; 
    
    lb(n+3)=0;%-ac; 
   
    ub(n+1)=15*pi;%init_guess(n+1)+4*pi;%0.26/30;%/30; 
    ub(n+4)=5;
    
    ub(n+2)=max_V; 
    ub(n+5)=max_V; 
    
    ub(n+3)=ac; 
    
end
%% Build the jacobian
t(1)=time_vec(1);
t(2)=tfin(2);
csi=(tfin-tin)/2;
A=[1 t(1) t(1)^2 t(1)^3; 1 t(2)^1 t(2)^2 t(2)^3; 0 1 2*t(1) 3*t(1)^2; 0 1 2*t(2) 3*t(2)^2 ];
A_inv=inv(A);
phi_d=[0 1 2*csi 3*csi^2]*A_inv;
phi=[1 csi csi^2 csi^3]*A_inv;
db_dxu=[1 0; 0 0; 0 0; 0 0];
df_dx=0;
du_dxu=[0 1];

d_residual=phi_d*db_dxu-(h/2)*(df_dx*phi*db_dxu+df_du*du_dxu);


%% solve
[x_sol,fval,eflag,fminoutput]=fmincon(cost_fun,init_guess,A,b,Aeq,Beq,lb,ub,nlcon,options);

%TODO impose constraint on thrust (magnitude<ac)
%%
% non_lin_con(init_guess,N,h,initial_states,final_states)
x_sol=init_guess;
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
plot(x_k,y_k)
figure(2)
bar(ux_k)

hold off

% function a = vecangle(v1,v2)
% n=[1,0,0]';
% x = cross(v1,v2);
% c = sign(dot(x,n)) * norm(x);
% a = atan2(c,dot(v1,v2));
% end