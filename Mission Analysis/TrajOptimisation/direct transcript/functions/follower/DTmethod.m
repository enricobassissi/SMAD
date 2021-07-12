function J = DTmethod(X, data)

%     N =data.n_int;
    %T = m;
    %time = data.time;
% alpha = m; beta = m;
    m(1)=X(7);
    m(2)=X(end-3);
%     
%     J = 0;
%     for k = 2:N-1
%         %integration of cost function
%         dt = time(k+1) - time(k);
%         J = J +  dt/6 * (T(k-1) + 4*T(k) + T(k + 1));
%     end
%     J = J + T(end)*dt;
%     T_inplane = a_in .*m * 1000;
%     T_outplane = a_out .*m * 1000;
%     T = (T_inplane.^2 + T_outplane.^2).^0.5;
    J = (m(1) - m(2))*data.MU;

end