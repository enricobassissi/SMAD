function [ISP, T, m_dot] = qinetiq_t5_future(P_in)
% Qinetiq T5 Future Datasheet
P = [300; 425; 560; 680];
Th = [10e-3; 15e-3; 20e-3; 25e-3];
Isp = [2400; 2700; 2850; 3010];

% Thrust fitting
fun_T = fit(P, Th, 'poly1');
T = feval(fun_T, P_in); % in mN

% Isp fitting
fun_ISP = fit(Th, Isp, 'power1'); % pretty close
ISP = feval(fun_ISP, T); % in s

% mass flow rate consequence
m_dot = T./ISP/9.8065; % in kg/s
end