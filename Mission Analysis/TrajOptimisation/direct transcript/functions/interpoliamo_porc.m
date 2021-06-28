function [SC2] = interpoliamo_porc(SC,n)

method = 'linear';
% method = 'spline';

t_new = linspace(SC.timead(1),SC.timead(end),n)';
SC2.timead       = t_new;
SC2.HS.X = interp1(SC.timead,SC.HS.X,t_new,method);
SC2.HS.r = interp1(SC.timead,SC.HS.r,t_new,method);
SC2.HS.beta = interp1(SC.timead,SC.HS.beta,t_new,method);
SC2.HS.alpha = interp1(SC.timead,SC.HS.alpha,t_new,method);
SC2.HS.T = interp1(SC.timead,SC.HS.T,t_new,method);

SC2.RES.T = interp1(SC.timead,SC.RES.T,t_new,method);
SC2.RES.T_outplane = interp1(SC.timead,SC.RES.T_outplane,t_new,method);
SC2.RES.T_inplane = interp1(SC.timead,SC.RES.T_inplane,t_new,method);
SC2.RES.gamma = interp1(SC.timead,SC.RES.gamma,t_new,method);
SC2.RES.el = interp1(SC.timead,SC.RES.el,t_new,method);

SC2.Xad = interp1(SC.timead,SC.Xad,t_new,method);
SC2.Xpropad = interp1(SC.timead,SC.Xpropad,t_new,method);
SC2.Xpropad_DT = interp1(SC.timead,SC.Xpropad_DT,t_new,method);

SC2.PROPULSION.DistanceSC_Sun_magnitude = interp1(SC.timead,SC.PROPULSION.DistanceSC_Sun_magnitude,t_new,method);

SC2.EPS.T_transf_orbit = interp1(SC.timead,SC.EPS.T_transf_orbit,t_new,method);
SC2.EPS.R_cartesian = interp1(SC.timead,SC.EPS.R_cartesian,t_new,method);
SC2.EPS.time = interp1(SC.timead,SC.EPS.time,t_new,method);

end

