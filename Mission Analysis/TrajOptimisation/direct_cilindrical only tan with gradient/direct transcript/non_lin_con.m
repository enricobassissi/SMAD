function [constraints, constraints_eq,GRAD,GRADeq ] = non_lin_con(y,N,time_vec,d_residual)
%EQUALITY CONSTRAINTS
%imposing differential equation through residual
residual=[];
d_residual_vec={};

for i=1:N-1
    n=(i-1)*5;
    x_k=[y(n+1) y(n+4) y(n+2) y(n+5)]'; %teta r vt vr 
    u_k=[y(n+3) 0]'; %ut ur
    n=n+5;
    x_k1=[y(n+1) y(n+4) y(n+2) y(n+5)]';
    u_k1=[y(n+3) 0]';
    h=time_vec(i+1)-time_vec(i);
    %Collocation point
    x_kc=0.5*(x_k+x_k1)+(h/8)*(dynamic(x_k, u_k) - dynamic(x_k1, u_k1));
%     x_dot_kc=-(3/(2*h))*(x_k+x_k1)+(h/8)*(x_dot_k-x_dot_k1);
    u_kc=(u_k+u_k1)/2;
    delta=x_k-x_k1+(h/6)*(dynamic(x_k,u_k)+4*dynamic(x_kc,u_kc)+dynamic(x_k1,u_k1));
    for j=1:4
        residual=[residual;delta(j)];
    end
    d_residual_mat = double(generate_jacobian(x_k,u_k(1),x_kc,u_kc(1),x_k1,u_k1(1),time_vec(i),time_vec(i+1),d_residual));
    d_res_1{i}=d_residual_mat(:,1:5);
    d_res_2{i}=d_residual_mat(:,6:10);

end
GRAD=[];
GRADeq=[blkdiag(d_res_1{:}) zeros((N-1)*4,5)]'+[zeros((N-1)*4,5) blkdiag(d_res_2{:})]';

%imposing boundaries condition
% boundaries=ones(2,1);
% n=0;
% vec=[y(n+1); y(n+2); y(n+3); y(n+4); y(n+5); y(n+6)];
% if isequal(vec,initial_states)
%     boundaries(1)=0;
% end
% n=N;
% vec=[y(n+1); y(n+2); y(n+3); y(n+4); y(n+5); y(n+6)];
% if isequal(vec,final_states)
%     boundaries(2)=0;
% end


%add here any other constraints if needed

%...

%construct the nonlineq constraint vector for fmincon
% constraints_eq=[residual(:);boundaries(:)]';
constraints_eq=[residual(:)]';
% if flag==1
%     constraints_eq=[constraints_eq(:); 1e15]';
% end
    

%GENERAL CONSTRAINTS
constraints=[];

%gradient

        
        

end
