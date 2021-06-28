clear all
%% 2020 VV
H=27.28;
min_dens=1000;%km/m3
max_dens=5000;%kg/m3
alow=0.3;
dlow= 10^( 3.1236 - 0.5*log10(alow) - 0.2*H);%km
Vlow=(4/3)*pi*((dlow*1000/2)^3) %m3
min_mass=Vlow*min_dens
aup=0.05;;
dup= 10^( 3.1236 - 0.5*log10(aup) - 0.2*H);
Vup=(4/3)*pi*((dup*1000/2)^3)
max_mass=Vup*max_dens
vesc_low=sqrt(6.67*10^(-11)*2*min_mass/(dlow/2))
vesc_up=sqrt(6.67*10^(-11)*2*max_mass/(dup/2))
%% 2009 TD 17
H=27.7;
min_dens=1000;%km/m3
max_dens=5000;%kg/m3
alow=0.3;
dlow= 10^( 3.1236 - 0.5*log10(alow) - 0.2*H);%km
Vlow=(4/3)*pi*((dlow*1000/2)^3) %m3
min_mass=Vlow*min_dens
aup=0.05;;
dup= 10^( 3.1236 - 0.5*log10(aup) - 0.2*H);
Vup=(4/3)*pi*((dup*1000/2)^3)
max_mass=Vup*max_dens
vesc_low=sqrt(6.67*10^(-11)*2*min_mass/(dlow/2))
vesc_up=sqrt(6.67*10^(-11)*2*max_mass/(dup/2))
%% 2011BP40
H=25.7;
min_dens=1000;%km/m3
max_dens=5000;%kg/m3
alow=0.3;
dlow= 10^( 3.1236 - 0.5*log10(alow) - 0.2*H);%km
Vlow=(4/3)*pi*((dlow*1000/2)^3) %m3
min_mass=Vlow*min_dens
aup=0.05;;
dup= 10^( 3.1236 - 0.5*log10(aup) - 0.2*H);
Vup=(4/3)*pi*((dup*1000/2)^3)
max_mass=Vup*max_dens
vesc_low=sqrt(6.67*10^(-11)*2*min_mass/(dlow/2))
vesc_up=sqrt(6.67*10^(-11)*2*max_mass/(dup/2))
%% 2021 JE 1
H=26.67;
min_dens=1000;%km/m3
max_dens=5000;%kg/m3
alow=0.3;
dlow= 10^( 3.1236 - 0.5*log10(alow) - 0.2*H);%km
Vlow=(4/3)*pi*((dlow*1000/2)^3) %m3
min_mass=Vlow*min_dens
aup=0.05;;
dup= 10^( 3.1236 - 0.5*log10(aup) - 0.2*H);
Vup=(4/3)*pi*((dup*1000/2)^3)
max_mass=Vup*max_dens
vesc_low=sqrt(6.67*10^(-11)*2*min_mass/(dlow/2))
vesc_up=sqrt(6.67*10^(-11)*2*max_mass/(dup/2))
%% grav acc
% G=6.674e-11; %m3/(kg*s2)
% mass= 7.59e9;%kg2.48e6;
% mu=mass*G;
% for i=1:10000
%     r(i)=200+i*3;
%     a(i)=mu/(r(i)^2);
% end
% plot(r,a)
% clear a
% %% circular orbits
% G=6.674e-11; %m3/(kg*s2)
% mass= 7.59e9;%kg2.48e6;
% mu=mass*G;
% for i=1:1000 %m
% a(i)=200+(i-1)*20;
% %let's try a circolar orbit with r=2km
% v(i)=sqrt(mu/a(i));%m/s
% T(i)=2*pi*sqrt(a(i)^3/mu); %s
% t_days(i)=T(i)/(3600*24);
% end
% plot(a,v)
% %% Perturbations
% S=1367;%W/m^2
% Cr=0.3;%coeff of reflectivity
% As=10;%m^2
% c=299792.458*1000%m/s
% Fsrp=S*Cr*As/c;
