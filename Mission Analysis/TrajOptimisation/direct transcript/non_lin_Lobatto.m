function [constraints, constraints_eq] = non_lin_Lobatto(y,N,time_vec)
%EQUALITY CONSTRAINTS
%imposing differential equation through residual
residual=[];
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
    
    x_1=(1/686)*((39*sqrt(21)+231)*x_k+224*x_kc+(-39*sqrt(21)+231)*x_k1+h*((3*sqrt(21)+21)*dynamic(x_k,u_k)-16*sqrt(21)*dynamic(x_kc,u_kc)+(3*sqrt(21)-21)*dynamic(x_k1,u_k1)));
    x_2=(1/686)*((-39*sqrt(21)+231)*x_k+224*x_kc+(39*sqrt(21)+231)*x_k1+h*((-3*sqrt(21)+21)*dynamic(x_k,u_k)+16*sqrt(21)*dynamic(x_kc,u_kc)+(-3*sqrt(21)-21)*dynamic(x_k1,u_k1)));
    u_1=u_k+(u_kc-u_k)*(1-sqrt(3/7));
    u_2=u_kc+(u_k1-u_kc)*sqrt(3/7);
    delta_1=(1/360)*(32*sqrt(21)+180)*x_k-64*sqrt(21)*x_kc+(32*sqrt(21)-180)*x_k1+h*((9+sqrt(21))*dynamic(x_k,u_k)+98*dynamic(x_1,u_1)+64*dynamic(x_kc,u_kc)+(9-sqrt(21))*dynamic(x_k1,u_k1));
    delta_2=(1/360)*(-32*sqrt(21)+180)*x_k+64*sqrt(21)*x_kc+(-32*sqrt(21)-180)*x_k1+h*((9-sqrt(21))*dynamic(x_k,u_k)+98*dynamic(x_2,u_2)+64*dynamic(x_kc,u_kc)+(9+sqrt(21))*dynamic(x_k1,u_k1));
    %delta=x_k-x_k1+(h/6)*(dynamic(x_k,u_k)+4*dynamic(x_kc,u_kc)+dynamic(x_k1,u_k1));
    for j=1:4
        residual=[residual;delta_1(j);delta_2(j)];
    end

end

%add here any other constraints if needed

%...

%construct the nonlineq constraint vector for fmincon
% constraints_eq=[residual(:);boundaries(:)]';
constraints_eq=[residual(:)]';

    

%GENERAL CONSTRAINTS
constraints=[];

end
