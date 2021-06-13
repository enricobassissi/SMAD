function [jacob_] = generate_jacobian(xt,ut,xct,uct,xt_1,ut_1,t_1,t_2,d_residual)
h1=(t_2-t_1);
jacob=d_residual(t_1,t_2,xt(1),xt(2),xt(3),xt(4),ut(1),xct(1),xct(2),xct(3),xct(4),uct(1), xt_1(1), xt_1(2), xt_1(3), xt_1(4), ut_1(1),h1);
jacob_=[jacob(:,1) jacob(:,3) jacob(:,5) jacob(:,2) jacob(:,4) jacob(:,6) jacob(:,8) jacob(:,10) jacob(:,7) jacob(:,9)] ;
end

