function M_ratio = thrust_integrator(y,N,t_vec,M0,tin,tfin,sim)
for i=1:N
    n=(i-1)*5;
%     u(i,:)=[y(n+7),y(n+8),y(n+9)];
    u(i)=abs(y(n+3));

end


acc=@(t) acc_stepwise(u,t,t_vec);
% for i=1:N
%     output.acc_vec(i)=acc(t_vec(i));
% end
[T,m]=ode45(@(t, mass) dynamic_mass(t,mass,acc,sim),linspace(0,tfin-tin,N) ,M0);
M_ratio=(m(1)-m(N))/m(1);


%%%
%         m = zeros(N,1); %% actually k - what is k? 
%         m(1) = M0;
%         K = -a_vect/(PS.Is*sim.g0)./l_d;
%         m(2) = m(1) +dl/2*(K(1)*m(1)+K(2)*m(2));
%         m(3) = m(2) +dl/2*(K(2)*m(2)+K(3)*m(3));
%         
%         for i = 3:sim.n_sol-1 % dal quarto al penultimo punto con predictor corrector AB3AM4 expl
%             m(i+1) = m(i) + dl12*(23*K(i)*m(i) - 16*K(i-1)*m(i-1) + 5*K(i-2)*m(i-2) ) ;
%             m(i+1) = m(i) + dl24*(9*K(i+1)*m(i+1)+19*K(i)*m(i) - 5*K(i-1)*m(i-1) + K(i-2)*m(i-2) ) ;
%         end
%%%
end

function K = dynamic_mass(t,mass,acc,sim) 
    acc_val=acc(t);
    Isp=sim.PS.Isp;
    g=sim.g0;
    K=-acc_val*mass/(Isp*g);
end

function out = acc_stepwise(u,t,t_vec)
%INPUT u=ux1 uy1 uz1 ux2 uy2 uz2 ...]
for i=1:(length(t_vec)-1)
    if and(t_vec(i)<=t,t<t_vec(i+1))
        break
    end
end
out=abs(u(i));
end