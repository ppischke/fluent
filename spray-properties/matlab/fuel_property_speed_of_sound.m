function [c_sound]= fuel_property_speed_of_sound(p,fuel)
T=300-273.15;
p=p/(10^5);
if strcmp('Butanol',fuel)
c_sound=1570.8;
elseif strcmp('Dodecane',fuel)
c_sound=1627.7;
elseif strcmp('Ethanol',fuel)
c_sound=1502.3;
elseif strcmp('iso-Octane',fuel)
c_sound=1502.6;
elseif strcmp('THFA',fuel)
c_sound=1706.1;
elseif strcmp('Heptane',fuel)
c_sound=1529.2;
elseif strcmp('Decanol',fuel)
c_sound=1682.7;
elseif strcmp('DNBE',fuel)
c_sound=NaN;
elseif strcmp('Ethyl-Levulinate',fuel)
    A1=1001;
    A2=-3.495;
    A3=-5.783*10^(-4);
    A4=19.02;
    A5=508;
    A6=4.576*10^(-4);
    c_sound=A1+A2*T+A3*p+A4*(p+A5)^(0.5)+A6*T*p; %m/s
elseif strcmp('Butyl-Levulinate',fuel)
    A1=450.2;
    A2=-3.19;
    A3=-0.1275;
    A4=32.99;
    A5=840.8;
    A6=4.346*10^(-4);
    c_sound=A1+A2*T+A3*p+A4*(p+A5)^(0.5)+A6*T*p; %m/s
elseif strcmp('2-MTHF',fuel)
    A1=420.7;
    A2=-3.712;
    A3=-0.1214;
    A4=33.12;
    A5=672.7;
    A6=5.827*10^(-4);
    c_sound=A1+A2*T+A3*p+A4*(p+A5)^(0.5)+A6*T*p; %m/s
elseif strcmp('DMS-T-1.5',fuel)
c_sound=1086.71;
elseif strcmp('DMS-T-2',fuel)
c_sound=1064.93;
elseif strcmp('DMS-T-5',fuel)
c_sound=1064.49;
elseif strcmp('DMS-T-10',fuel)
c_sound=1066.34;
elseif strcmp('Diesel',fuel)
c_sound=NaN;
end
end