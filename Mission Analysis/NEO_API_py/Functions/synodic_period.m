function S = synodic_period(a1, a2)

    muSun = astroConstants(4);
    T1 = 2*pi*sqrt(a1^3/muSun);
    T2 = 2*pi*sqrt(a2^3/muSun);
    
    if T1 < T2
        S = 1/(1/T1 - 1/T2);
    else
        S = 1/(1/T2 - 1/T1);
    end

end