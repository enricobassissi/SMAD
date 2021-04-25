function [dvtot,VI,VF,ToF]=lambert_solver_flyby2( R1, R2, V1, V2, t1, t2, ksun, vlim)
%{
CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John
%}
ToF = t2 - t1;
[a,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(R1,R2,ToF,ksun,0,0,0,0);

dv=sqrt((VF(1)-V2(1))^2+(VF(2)-V2(2))^2+(VF(3)-V2(3))^2);

dv=abs(dv);

if dv>vlim
    dvtot=NaN;
else
    dvtot=dv;
end

end