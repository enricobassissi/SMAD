function [dvtot,VI,VF,ToF]=lambert_solver_flyby( R1, R2, V1, V2, t1, t2, ksun, vlim)

%{
CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John
%}
ToF = t2 - t1;
[a,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(R1,R2,ToF,ksun,0,0,0,0);

dv=sqrt((VI(1)-V1(1))^2+(VI(2)-V1(2))^2+(VI(3)-V1(3))^2);

dv=abs(dv);

if dv>vlim
    dvtot=NaN;
else
    dvtot=dv;
end

end