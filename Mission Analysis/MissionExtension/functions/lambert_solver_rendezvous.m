function [dvtot,dv1,dv2]=lambert_solver_rendezvous( R1, R2, V1, V2, TOF, ksun)
% tof in sec
[~,~,~,~,VI,VF,~,~] = lambertMR(R1,R2,TOF,ksun,0,0,0,0);

dv1=sqrt((VI(1)-V1(1))^2+(VI(2)-V1(2))^2+(VI(3)-V1(3))^2);
dv2=sqrt((VF(1)-V2(1))^2+(VF(2)-V2(2))^2+(VF(3)-V2(3))^2);

dvtot = abs(dv1)+abs(dv2);


end