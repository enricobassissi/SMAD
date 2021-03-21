function S = synodic_period(pl1,pl2)
%{
CONTRIBUTORS: Apparenza Lucia
              Bassissi Enrico
              Di Trocchio Marco
              Lane John
%}
date_A =  [2020, 1, 1, 0, 0, 0];
T = date2mjd2000(date_A);
            
switch pl1
    case 'mercury'
        if pl2=='mars'
            [kep_MERC,ksun] = uplanet(T, 1);
            P_MER=2*pi*sqrt(kep_MERC(1)^3/ksun);
            
            [kep_MARS,ksun] = uplanet(T, 4);
            P_MAR=2*pi*sqrt(kep_MARS(1)^3/ksun);
            
            S=(P_MAR*P_MER)/abs((P_MER-P_MAR));
        end
    case 'mars'
        if pl2=='neptune'
            [kep_NEP,ksun] = uplanet(T, 8);
            P_NEP=2*pi*sqrt(kep_NEP(1)^3/ksun);
            
            [kep_MARS,ksun] = uplanet(T, 4);
            P_MAR=2*pi*sqrt(kep_MARS(1)^3/ksun);
            
            S=(P_MAR*P_NEP)/abs((P_NEP-P_MAR));
        end
end
end