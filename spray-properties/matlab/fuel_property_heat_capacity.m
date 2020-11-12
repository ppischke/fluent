function [CP]= fuel_property_heat_capacity(T,fuel)

%if fuel==char('Butanol')
if strcmp('Butanol',fuel)
    A=1.912*10^5;
    B=-7.304*10^2;
    C=2.2998*10^0;
    D=0;
    M_fuel=74.12;
    CP = (A+B*T+C*T^2)/M_fuel;
elseif strcmp('Dodecane',fuel)
    A=5.0821*10^5;
    B=-1.3621*10^3;
    C=3.1015*10^0;
    D=0;
    M_fuel=170.34;
    CP = (A+B*T+C*T^2+D*T^3)/M_fuel;
elseif strcmp('Ethanol',fuel)
    A=1.0264*10^5;
    B=-1.3963*10^2;
    C=-3.0341*10^-2;
    D=2.0386*10^-3;
    M_fuel=46.07;
    CP = (A+B*T+C*T^2+D*T^3)/M_fuel;
elseif strcmp('iso-Octane',fuel)
    A=9.5275*10^4;
    B=6.967*10^2;
    C=-1.3765*10^0;
    D=2.1734*10^-3;
    M_fuel=114.23;
    CP = (A+B*T+C*T^2+D*T^3)/M_fuel;
elseif strcmp('THFA',fuel)
    A=5.27*10^4;
    B=4.358*10^2;
    C=0;
    D=0;
    M_fuel=102.13;
    CP = (A+B*T+C*T^2+D*T^3)/M_fuel;
elseif strcmp('Heptane',fuel)
    A=61.26;
    B=314410;
    C=1824.6;
    D=-2547.9;
    M_fuel=100.20;
    T_c=540.2;
    t=(1-T/T_c);
    CP = (A^2/t+B-2*A*C*t-A*D*t^2-C^2*t^3/3-C*D*t^4/2-D^2*t^5/5)/M_fuel;
elseif strcmp('Decanol',fuel)
    A=4988500;
    B=-52898;
    C=216.35;
    D=-0.37538;
    E=0.00023674;
    M_fuel=158.28;
    CP = (A+B*T+C*T^2+D*T^3+E*T^4)/M_fuel;
elseif strcmp('DNBE',fuel)
    A=2.7072E+5;
    B=-2.5983E+2;
    C=9.5427E-1;
    D=0;
    E=0;
    M_fuel=130.23;
    CP = (A+B*T+C*T^2+D*T^3+E*T^4)/M_fuel;
elseif strcmp('2-MTHF',fuel) %Kai Leonhardt Excel File
    A=5.1045E-04;
    B=5.6431E-03;
    C=1.0570E+02;
    M_fuel=86.1323;%g/mol
    CP = (A*T^2+B*T+C)/(M_fuel/1000);
end
end




