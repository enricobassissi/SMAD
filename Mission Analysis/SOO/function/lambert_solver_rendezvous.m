function [dvtot,dv1,dv2]=lambert_solver_rendezvous( R1, R2, V1, V2, t1, t2, ksun, vlim)
%[dvtot,VI,VF,TOF]=lambert_solver_rendezvous( R1, R2, V1, V2, t1, t2, ksun, vlim)
ToF_sec = t2 - t1;
% TOF = ToF_sec/(3600*24);
[~,~,~,~,VI,VF,~,~] = lambertMR(R1,R2,ToF_sec,ksun,0,0,0,0);

dv1=sqrt((VI(1)-V1(1))^2+(VI(2)-V1(2))^2+(VI(3)-V1(3))^2);
dv2=sqrt((VF(1)-V2(1))^2+(VF(2)-V2(2))^2+(VF(3)-V2(3))^2);

dv = abs(dv1)+abs(dv2);

if dv>vlim
    dvtot=NaN;
else
    dvtot=dv;
end

end