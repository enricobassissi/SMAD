function [] = DT_plots(HS, RES, timead, Xad, Xpropad, Xpropad_DT, IDcase, colors)
% --- plot check of the variables
% r, TH, z, vr, theta_dot, vz, m
% figure('Name','actual state vs Propagated')
% sgtitle(string({'Comparison NLI - State Propagation - ',IDcase}))
% subplot(7,1,1)
% plot(timead, Xad(:,1),'Color', colors(1,:));
% hold on; grid on; ylabel('r [AU]');
% plot(timead, Xpropad(:,1),'Color', colors(2,:));
% subplot(7,1,2)
% plot(timead, Xad(:,2),'Color', colors(1,:));
% hold on; grid on; ylabel('$\theta$ [rad]');
% plot(timead, Xpropad(:,2),'Color', colors(2,:));
% subplot(7,1,3)
% plot(timead, Xad(:,3),'Color', colors(1,:));
% hold on; grid on; ylabel('z [AU]');
% plot(timead, Xpropad(:,3),'Color', colors(2,:));
% subplot(7,1,4)
% plot(timead, Xad(:,4),'Color', colors(1,:));
% hold on; grid on; ylabel('$v_r$ [-]');
% plot(timead, Xpropad(:,4),'Color', colors(2,:));
% subplot(7,1,5)
% plot(timead,Xad(:,5),'Color', colors(1,:));
% hold on; grid on; ylabel('$\dot{\theta}$ [-]');
% plot(timead,Xpropad(:,5),'Color', colors(2,:));
% subplot(7,1,6)
% plot(timead,Xad(:,6),'Color', colors(1,:));
% hold on; grid on; ylabel('$v_z$ [-]');
% plot(timead,Xpropad(:,6),'Color', colors(2,:));
% subplot(7,1,7)
% plot(timead,Xad(:,7),'Color', colors(1,:));
% hold on; grid on; ylabel('$m_f$ [-]');
% plot(timead,Xpropad(:,7),'Color', colors(2,:))
% legend('NLI','Prop')

% % ----- Adimensionalization of variables
% figure('Name','Trajectory parameters'),
% subplot(5,1,1), plot(timead,rad2deg(RES.el)), title('out of plane thrust angle'), grid on
% subplot(5,1,2), plot(timead,rad2deg(RES.gamma)), title('in plane thrust angle'), grid on
% subplot(5,1,3), plot(timead,RES.T_inplane), title('T inplane'), grid on
% subplot(5,1,4), plot(timead,RES.T_outplane), title('T outplane'), grid on
% subplot(5,1,5), plot(timead,RES.T), title('T'), grid on

% ----- control allocation
figure('Name','Control Allocation')
sgtitle(string({'Control Allocation - ',IDcase}))
subplot(3,1,1)
plot(timead, RES.T,'Color', colors(1,:))
hold on; grid on; title('T [N]')
plot(timead, HS.T,'Color', colors(2,:))
subplot(3,1,2)
plot(timead,rad2deg(RES.gamma),'Color', colors(1,:))
hold on, grid on; title('in plane thrust angle [deg]')
plot(timead, rad2deg(HS.alpha),'Color', colors(2,:))
subplot(3,1,3)
plot(timead, rad2deg(RES.el),'Color', colors(1,:))
hold on; grid on; title('out of plane thrust angle [deg]')
plot(timead, rad2deg(wrapToPi(HS.beta)),'Color', colors(2,:))
legend('NLI','DT')

% ----- CHECK PHYSICS VALIDITY
figure('Name','NLI vs Propagated DT')
sgtitle(string({'Comparison NLI - DT Propagation - ',IDcase}))
subplot(7,1,1)
plot(timead, Xad(:,1),'Color', colors(1,:))
hold on; grid on; ylabel('r [AU]');
plot(timead, Xpropad_DT(:,1),'Color', colors(2,:))
subplot(7,1,2)
plot(timead, Xad(:,2),'Color', colors(1,:))
hold on; grid on; ylabel('$\theta$ [rad]');
plot(timead, Xpropad_DT(:,2),'Color', colors(2,:))
subplot(7,1,3)
plot(timead, Xad(:,3),'Color', colors(1,:))
hold on; grid on; ylabel('z [AU]');
plot(timead, Xpropad_DT(:,3),'Color', colors(2,:))
subplot(7,1,4)
plot(timead, Xad(:,4),'Color', colors(1,:))
hold on; grid on; ylabel('$v_r$ [-]');
plot(timead, Xpropad_DT(:,4),'Color', colors(2,:))
subplot(7,1,5)
plot(timead, Xad(:,5),'Color', colors(1,:))
hold on; grid on; ylabel('$\dot{\theta}$ [-]');
plot(timead, Xpropad_DT(:,5),'Color', colors(2,:))
subplot(7,1,6)
plot(timead, Xad(:,6),'Color', colors(1,:))
hold on; grid on; ylabel('$v_z$ [-]');
plot(timead, Xpropad_DT(:,6),'Color', colors(2,:))
subplot(7,1,7)
plot(timead, Xad(:,7),'Color', colors(1,:))
hold on; grid on; ylabel('$m_f$ [-]');
plot(timead, Xpropad_DT(:,7),'Color', colors(2,:))
legend('NLI','DT')