function output = interpoliamo_insieme(output)  

t_new = linspace(output.t(1),output.t(end),length(output.t));

output.u_equi = interp1(output.t,output.u,t_new);
output.m_equi = interp1(output.t,output.m,t_new);
output.s_equi = interp1(output.t,output.s,t_new);
output.r_equi = interp1(output.t,output.r,t_new);
output.l_equi = interp1(output.t,output.l,t_new);
output.z_equi = interp1(output.t,output.z,t_new);
output.v_r_equi = interp1(output.t,output.v_r,t_new);
output.w_equi = interp1(output.t,output.w,t_new);
output.v_z_equi = interp1(output.t,output.v_z,t_new);

end


