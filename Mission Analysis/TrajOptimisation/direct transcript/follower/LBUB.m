function [lb, ub] = LBUB(XX0, X,  data)

N = data.n_int;
T_lb = 0;
T_ub = data.Tmax;
alpha_bound = data.angle_inplane_panels_max;
% gamma, in plane thrust angle
alpha_lb = deg2rad(-alpha_bound);
alpha_ub = deg2rad(alpha_bound);
% out of plane thrust angle
beta_lb = deg2rad(-360); % JP style
beta_ub = deg2rad(360);

%limits definition
ub = zeros(1,10*N);
lb = zeros(1,10*N);

% -- dynamic boundaries : sigh not now
%     XX0_dyn_bound_ub = reshape(XX0,10,50); % a little bigger than XX0
%     XX0_dyn_bound_lb = reshape(XX0,10,50); % a little smaller than XX0

%r
r_ub = 2;
r_lb = 0.7;
%     r_ub = max(abs(XX0_dyn_bound_ub(1,:)))*1.5;%     
%     r_lb = min(abs(XX0_dyn_bound_lb(1,:)))*0.5;%     
%theta
th_ub = 20;
th_lb = 0;
%     th_ub = max(XX0_dyn_bound_ub(2,:)); %     
%     th_lb = min(XX0_dyn_bound_lb(2,:));%     
%z
z_ub =  0.05;
z_lb =  -0.05;
% 	z_ub = max(XX0_dyn_bound_ub(3,:));%     
%     z_lb = sign(abs(min(XX0_dyn_bound_lb(3,:)))*1.5);%     

%vr
vr_ub = 0.3;
vr_lb = -0.3;
%     vr_ub = max(XX0_dyn_bound_ub(4,:));% 
%     vr_lb = min(XX0_dyn_bound_lb(4,:));% 
%theta_dot
thd_ub = 5;
thd_lb = 0;
%     thd_ub = max(XX0_dyn_bound_ub(5,:));%
%     thd_lb = min(XX0_dyn_bound_ub(5,:));%
%vz
vz_ub = 0.03;
vz_lb = -0.03;
%     vz_ub = max(XX0_dyn_bound_ub(6,:));%
%     vz_lb = min(XX0_dyn_bound_lb(6,:));%

%m
m_ub = X(1,7);
m_lb = 0.7;
%     m_lb = min(XX0_dyn_bound_lb(7,:));

%INITIAL CONDITION
ub(1:10) = XX0(1:10); 
ub(8) = T_ub; ub(9) = alpha_ub; ub(10) = beta_ub;
ub(7) = m_ub;
lb(1:10) = XX0(1:10); 
lb(7) = m_ub; % the initial mass is fixed, so its lower/upper boundaries are the same
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
lb(17:10:((N-1)*10)) = m_lb;
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