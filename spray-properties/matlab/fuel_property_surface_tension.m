function [surface_tension]= fuel_property_surface_tension(T,fuel)

if strcmp('Butanol',fuel)
    A=4.983*10^-2;
    B=-8.54*10^-5;
    C=0;
    D=0;
    E=0;
    surface_tension = A+B*T+C*T^2+D*T^3+E*T^4;
elseif strcmp('Dodecane',fuel)
    A=5.5493*10^-2;
    B=1.3262*10^0;
    C=0;
    D=0;
    E=0;
    T_r=T/658;
    surface_tension = A*(1-T_r)^(B+C*T_r+D*T_r^2+E*T_r^3);
elseif strcmp('Ethanol',fuel)
    A=3.764*10^-2;
    B=-2.157*10^-5;
    C=-1.025*10^-7;
    D=0;
    E=0;
    surface_tension = A+B*T+C*T^2+D*T^3+E*T^4;
elseif strcmp('iso-Octane',fuel)
    A=4.7675*10^-2;
    B=1.2018*10^0;
    C=0;
    D=0;
    E=0;
    T_r=T/543.8;
    surface_tension = A*(1-T_r)^(B+C*T_r+D*T_r^2+E*T_r^3);
elseif strcmp('THFA',fuel)
    A=6.5073*10^-2;
    B=-8.8258*10^-5;
    C=-1.6138*10^-8;
    D=0;
    E=0;
    surface_tension = A+B*T+C*T^2+D*T^3+E*T^4;
elseif strcmp('Heptane',fuel)
    A=0.054143;
    B=1.2512;
    C=0;
    D=0;
    E=0;
    T_r=T/540.2;
    surface_tension = A*(1-T_r)^(B+C*T_r+D*T_r^2+E*T_r^3);
elseif strcmp('Decanol',fuel)
    A=0.051263;
    B=1.0395;
    C=0;
    D=0;
    E=0;
    T_r=T/687.3;
    surface_tension = A*(1-T_r)^(B+C*T_r+D*T_r^2+E*T_r^3);
elseif strcmp('DNBE',fuel)
    A=5.1346E-02;
    B=1.1604E+00;
    C=0;
    D=0;
    E=0;
    T_r=T/584.1;
    surface_tension = A*(1-T_r)^(B+C*T_r+D*T_r^2+E*T_r^3);
elseif strcmp('DMS-T-1.5',fuel)
    surface_tension = 17.31039/1000;
elseif strcmp('DMS-T-2',fuel)
    surface_tension = 18.2586/1000;
elseif strcmp('DMS-T-5',fuel)
    surface_tension = 18.3358/1000;
elseif strcmp('DMS-T-10',fuel)
    surface_tension = 19.51667/1000;
elseif strcmp('Ethyl-Levulinate',fuel)
    A=-0.1045;
    B=63.87;
    surface_tension = (A*T+B)/1000; %http://pubs.rsc.org/en/content/articlehtml/2011/gc/c0gc00853b
elseif strcmp('Butyl-Levulinate',fuel)
    A=-0.0902;
    B=57.83;
    surface_tension = (A*T+B)/1000; %http://pubs.rsc.org/en/content/articlehtml/2011/gc/c0gc00853b
elseif strcmp('2-MTHF',fuel)
    surface_tension = 24.6/1000; %J. of Chem. Eng Data Vol. 50 S. 1334-1337
elseif strcmp('Diesel',fuel)
    surface_tension = 22.7/1000;
end
end




