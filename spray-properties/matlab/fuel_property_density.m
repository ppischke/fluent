function [density]        = fuel_property_density(T,fuel);

if strcmp('Butanol',fuel)
    A=9.65*10^-1;
    B=2.666*10^-1;
    C=5.6305*10^2;
    D=2.4419*10^-1;
    M_fuel=74.12;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
    
elseif strcmp('Dodecane',fuel)
    A=3.5541*10^-1;
    B=2.5511*10^-1;
    C=6.58*10^2;
    D=2.9368*10^-1;
    M_fuel=170.34;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
    
elseif strcmp('Ethanol',fuel)
    A=1.648*10^0;
    B=2.7627*10^-1;
    C=5.1392*10^2;
    D=2.331*10^-1;
    M_fuel=46.07;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
    
elseif strcmp('iso-Octane',fuel)
    A=5.9059*10^-1;
    B=2.7424*10^-1;
    C=5.438*10^2;
    D=2.847*10^-1;
    M_fuel=114.23;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
    
elseif strcmp('THFA',fuel)
    A=9.7001*10^-1;
    B=2.813*10^-1;
    C=6.39*10^2;
    D=2.3837*10^-1;
    M_fuel=102.13;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
    
elseif strcmp('Heptane',fuel)
    A=0.61259;
    B=0.26211;
    C=540.2;
    D=0.28141;
    M_fuel=100.20;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
    
elseif strcmp('Decanol',fuel)
    A=0.37384;
    B=0.24241;
    C=687.3;
    D=0.26646;
    M_fuel=158.28;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;
   
elseif strcmp('Ethyl-Levulinate',fuel)
    A=922.4;
    B=-0.764;
    C=-0.01008;
    D=3.866;
    E=695;
    F=1.509*10^-4;
    p=1;
    TC=T-273.15;
    density = A + B*TC+C*p+D*(p+E).^(0.5)+F*TC*p;
    
elseif strcmp('Butyl-Levulinate',fuel)
    A=856.2;
    B=-0.7779;
    C=-0.02112;
    D=4.48;
    E=724.4;
    F=1.346*10^-4;
    p=1;
    TC=T-273.15;
    density = A + B*TC+C*p+D*(p+E).^(0.5)+F*TC*p;
    
elseif strcmp('2-MTHF',fuel)
    A=703.9;
    B=-0.823;
    C=-0.02625;
    D=5.778;
    E=762.4;
    F=1.559*10^-4;
    p=1;
    TC=T-273.15;
    density = A + B*TC+C*p+D*(p+E).^(0.5)+F*TC*p;
    
elseif strcmp('DMS-T-1.5',fuel)
    density = 848.99;
    
elseif strcmp('DMS-T-2',fuel)
    density = 887.087;

elseif strcmp('DMS-T-5',fuel)
    density = 905.99;

elseif strcmp('DMS-T-10',fuel)
    density = 933;

elseif strcmp('DNBE',fuel)
    A=5.5941E-01;
    B=2.7243E-01;
    C=5.8410E+02;
    D=2.9932E-01;
    M_fuel=130.23;
    density = A/(B^(1+(1-T/C)^D))*M_fuel;

elseif strcmp('Diesel',fuel)
    density = 833;

end
end




