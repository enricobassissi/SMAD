% --- times
SC1.uniform.time = [SC1.mjd2000_dep+SC1.leg1.timead.*sim.TU/86400; SC1.coasting.leg1.time;
    SC1.coasting.leg1.time(end)+SC1.leg2.timead.*sim.TU/86400; SC1.coasting.leg2.time;]; % mjd 2000
SC2.uniform.time = [SC2.mjd2000_dep+SC2.lega.timead.*sim.TU/86400; SC2.coasting.lega.time;
    SC2.coasting.lega.time(end)+SC2.legb.timead.*sim.TU/86400; SC2.coasting.legb.time;]; % mjd 2000

time_vect_SC1 = SC1.uniform.time;
time_vect_SC1(100) = (time_vect_SC1(101)+time_vect_SC1(99))/2;
time_vect_SC1(201) = (time_vect_SC1(202)+time_vect_SC1(200))/2;
time_vect_SC1(299) = time_vect_SC1(298)+(time_vect_SC1(301)-time_vect_SC1(298))/3;
time_vect_SC1(300) = time_vect_SC1(298)+(time_vect_SC1(301)-time_vect_SC1(298))*2/3;

time_vect_SC2 = SC2.uniform.time;
time_vect_SC2(99) = time_vect_SC2(98)+(time_vect_SC2(101)-time_vect_SC2(98))/3;
time_vect_SC2(100) = time_vect_SC2(98)+(time_vect_SC2(101)-time_vect_SC2(98))*2/3;
time_vect_SC2(201) = (time_vect_SC2(202)+time_vect_SC2(200))/2;
time_vect_SC2(300) = (time_vect_SC2(301)+time_vect_SC2(299))/2;

for i=1:length(time_vect_SC1)
    time_vect_plot_SC1(i) = datenum(datetime(mjd20002date(time_vect_SC1(i))));
    time_vect_plot_SC2(i) = datenum(datetime(mjd20002date(time_vect_SC2(i))));
end

% -- requested thrust
SC1.uniform.Thrust = [SC1.leg1.EPS.T_transf_orbit; zeros(data.n_int,3);
    SC1.leg2.EPS.T_transf_orbit; zeros(data.n_int,3)]; % N
SC2.uniform.Thrust = [SC2.lega.EPS.T_transf_orbit; zeros(data.n_int,3);
    SC2.legb.EPS.T_transf_orbit; zeros(data.n_int,3)]; % N

% -- available thrust
SC1.leg1.PROPULSION = nikita(SC1.leg1, data, sim, colors);
SC1.leg2.PROPULSION = nikita(SC1.leg2, data, sim, colors);
SC2.lega.PROPULSION = nikita(SC2.lega, data, sim, colors);
SC2.legb.PROPULSION = nikita(SC2.legb, data, sim, colors);
SC1.available.Thrust = [SC1.leg1.PROPULSION.T_max_modulated; zeros(data.n_int,1);
    SC1.leg2.PROPULSION.T_max_modulated; zeros(data.n_int,1)]; % N
SC2.available.Thrust = [SC2.lega.PROPULSION.T_max_modulated; zeros(data.n_int,1);
    SC2.legb.PROPULSION.T_max_modulated; zeros(data.n_int,1)]; % N

% -- NLI thrust
SC1.NLI.Thrust = [SC1.leg1.RES.T; zeros(data.n_int,1); SC1.leg2.RES.T; zeros(data.n_int,1)]; % N
SC2.NLI.Thrust = [SC2.lega.RES.T; zeros(data.n_int,1); SC2.legb.RES.T; zeros(data.n_int,1)]; % N

% -- plot
idx_sel = [1,25,50,75,100,225,250,275,400];
figure('Name','T SC1')
subplot(2,1,1)
plot(time_vect_plot_SC1, SC1.NLI.Thrust*1e3,'Color', colors(3,:),'DisplayName','T NLI')
hold on; grid on;
plot(time_vect_plot_SC1, SC1.available.Thrust*1e3,'--','Color', colors(2,:),...
    'DisplayName','T available','LineWidth',2.5)
plot(time_vect_plot_SC1, vecnorm(SC1.uniform.Thrust,2,2)*1e3,'Color', colors(1,:),'DisplayName','T requested')
legend('show')
% xlabel('Date Dep'); 
ylabel('T magnitude [mN]');
ax = gca;
ax.XTick=time_vect_plot_SC1(idx_sel) ;
xtickangle(30)
datetick('x','yyyy mmm dd','keepticks')
xlim ([time_vect_plot_SC1(1) time_vect_plot_SC1(end)])

subplot(2,1,2)
plot(time_vect_plot_SC2, SC2.NLI.Thrust*1e3,'Color', colors(3,:),'DisplayName','T NLI')
hold on; grid on;
plot(time_vect_plot_SC2, SC2.available.Thrust*1e3,'--','Color', colors(2,:),...
    'DisplayName','T available','LineWidth',2.5)
plot(time_vect_plot_SC2, vecnorm(SC2.uniform.Thrust,2,2)*1e3,'Color', colors(1,:),'DisplayName','T requested')
legend('show')
xlabel('Date Dep'); ylabel('T magnitude [mN]');
ax = gca;
ax.XTick=time_vect_plot_SC2(idx_sel) ;
xtickangle(30)
datetick('x','yyyy mmm dd','keepticks')
xlim ([time_vect_plot_SC2(1) time_vect_plot_SC2(end)])