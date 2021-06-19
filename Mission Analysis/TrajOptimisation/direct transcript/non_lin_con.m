function [constraints, constraints_eq] = non_lin_con(y,N,time_vec)
%EQUALITY CONSTRAINTS
%imposing differential equation through residual
residual=[];
flag=0;
for i=1:N-1
    n=(i-1)*5;
    x_k=[y(n+1) y(n+4) y(n+2) y(n+5)]';
    u_k=[y(n+3) 0]';
    n=n+5;
    x_k1=[y(n+1) y(n+4) y(n+2) y(n+5)]';
    u_k1=[y(n+3) 0]';
    h=time_vec(i+1)-time_vec(i);

    x_kc=0.5*(x_k+x_k1)+(h/8)*(dynamic(x_k, u_k) - dynamic(x_k1, u_k1));
%     x_dot_kc=-(3/(2*h))*(x_k+x_k1)+(h/8)*(x_dot_k-x_dot_k1);
    u_kc=(u_k+u_k1)/2;
    delta=x_k-x_k1+(h/6)*(dynamic(x_k,u_k)+4*dynamic(x_kc,u_kc)+dynamic(x_k1,u_k1));
    for j=1:4
        residual=[residual;delta(j)];
    end
%     if or((norm(u_k)>ac),norm(u_k1)>ac)
%         flag=1;
%     end
end
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

end
