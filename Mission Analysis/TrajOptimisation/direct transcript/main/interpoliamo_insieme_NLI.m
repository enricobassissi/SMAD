function output2 = interpoliamo_insieme_NLI(output,n)  

%{
% example:
t_old = linspace(0,1,10);
t_query = linspace(0,1,100)';
y_old = exp(2)*sin(3*t_old);
y_query = interp1(t_old,y_old,t_query,'linear');
plot(t_old,y_old,'o','DisplayName','Old'); 
hold on; 
plot(t_query,y_query,'.','DisplayName','New');
legend('show')

%%
% y = interp(x,r) increases the sample rate of x, the input signal, by a factor of r.
t = 0:1/1e3:1;
x = sin(2*pi*30*t) + sin(2*pi*60*t);
y = interp(x,4);
subplot(2,1,1)
stem(0:30,x(1:31),'filled','MarkerSize',3)
grid on
xlabel('Sample Number')
ylabel('Original')

subplot(2,1,2)
stem(0:120,y(1:121),'filled','MarkerSize',3)
grid on
xlabel('Sample Number')
ylabel('Interpolated')
%}

method = 'linear';

t_new = linspace(output.t(1),output.t(end),n)';
output2.t       = t_new;
output2.m       = interp1(output.t,output.m,t_new,method);
output2.Thrust  = interp1(output.t,output.Thrust,t_new,method);
output2.T_magn  = sqrt(output2.Thrust(:,1).^2 + output2.Thrust(:,3).^2);
output2.r       = interp1(output.t,output.r,t_new,method);
output2.theta   = interp1(output.t,output.theta,t_new,method);
output2.z       = interp1(output.t,output.z,t_new,method);
output2.a       = output.a;
output2.Href    = output.Href;
output2.vr      = interp1(output.t,output.vr,t_new,method);
output2.vt      = interp1(output.t,output.vt,t_new,method);
output2.vz      = interp1(output.t,output.vz,t_new,method);
output2.acc_inplane = interp1(output.t,output.acc_inplane,t_new,method);
output2.acc_out = interp1(output.t,output.acc_out,t_new,method);
output2.acc = interp1(output.t,output.acc,t_new,method);
output2.L = output.L;
output2.gamma = interp1(output.t,output.gamma,t_new,method);
output2.T_inplane = interp1(output.t,output.T_inplane,t_new,method);
output2.T_outplane = interp1(output.t,output.T_outplane,t_new,method);
output2.theta_dot = interp1(output.t,output.theta_dot,t_new,method);

end


