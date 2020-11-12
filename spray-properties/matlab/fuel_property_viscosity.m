function [viscosity]= fuel_property_viscosity(T,fuel)
if strcmp('Butanol',fuel)
    A=-35.426;
    B=3184.5;
    C=3.2965;
    D=-3E-27;
    E=10;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('Dodecane',fuel)
    A=-20.607;
    B=1943;
    C=1.3205;
    D=0;
    E=0;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('Ethanol',fuel)
    A=7.875;
    B=781.98;
    C=-3.0418;
    D=0;
    E=0;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('iso-Octane',fuel)
    A=-12.928;
    B=1137.5;
    C=0.25725;
    D=-3.6929E-28;
    E=10;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('THFA',fuel)
    A=-7.9742;
    B=2745.4;
    C=-1.1468;
    D=0;
    E=0;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('Heptane',fuel)
    A=-2.45E+01;
    B=1.53E+03;
    C=2.0087;
    D=0;
    E=0;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('Decanol',fuel)
    A=-8.07E+01;
    B=6.33E+03;
    C=9.646;
    D=0;
    E=0;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('DNBE',fuel)
    A=1.0027E+01;
    B=2.0600E+02;
    C=-3.1607E+00;
    D=0;
    E=0;
    viscosity = exp(A+B/T+C*log(T)+D*T^E);
elseif strcmp('Ethyl-Levulinate',fuel)
    A=-14.055;
    B=4.354*10^-4;
    C=12194;
    D=803.08;
    E=7.248*10^-3;
    p=1;
    TC=T-273.15;
    viscosity = exp(A+B*p+C/(TC+D)+E*(p/TC))/1000;
elseif strcmp('Butyl-Levulinate',fuel)
    A=-14.24;
    B=5.751*10^-4;
    C=12196;
    D=779.37;
    E=4.903*10^-3;
    p=1;
    TC=T-273.15;
    viscosity = exp(A+B*p+C/(TC+D)+E*(p/TC))/1000;
elseif strcmp('2-MTHF',fuel)
    A=-8.4912;
    B=3.227*10^-4;
    C=12169;
    D=1578.2;
    E=3.3034*10^-3;
    p=1;
    TC=T-273.15;
    viscosity = exp(A+B*p+C/(TC+D)+E*(p/TC))/1000;
elseif strcmp('DMS-T-1.5',fuel)
    viscosity = 1.25565/1000;
elseif strcmp('DMS-T-2',fuel)
    viscosity = 2.50646/1000;
elseif strcmp('DMS-T-5',fuel)
    viscosity = 4.47014/1000;
elseif strcmp('DMS-T-10',fuel)
    viscosity = 9.78095/1000;
elseif strcmp('Diesel',fuel)
    viscosity = 3.0/1000;
end
end