N=500;

x_= [x_sol(1),x_sol(4),x_sol(2),x_sol(5)]';%
% x_= [init_guess(1),init_guess(4),init_guess(2),init_guess(5)];%,sqrt(1/init_guess(4)),0];%
% 
clear u
u=ux_k;

acceleration= @(t) acc_stepwise(u,t,time_vec);
opts = odeset('Reltol',1e-12,'AbsTol',1e-14,'Stats','on');
[T,state_sol]=ode113(@(t,X) dynamic_checker(t,X, acceleration), linspace(tin, tfin,1000000), x_,opts);
%plot3(x_(:,1),x_(:,2),x_(:,3))
clear x_k; clear y_k
for i=1:(length(state_sol))
    x_k(i)=state_sol(i,1); 
    y_k(i)=state_sol(i,2); 
end
x_k=unwrap(x_k);
for i=1:(length(state_sol))
    [x_k_pol(i),y_k_pol(i)]=pol2cart(x_k(i),y_k(i));
end
figure(1)
plot(0,0,'x');
hold on
plot(x_k_pol,y_k_pol)
axis equal
hold off

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