function [B2]        = fuel_property_virial_coefficient(T,fuel);
%m3/kmol
if strcmp('Butanol',fuel)
    A=1.89E-1;
    B=-1.82E+2;
    C=-4.05E+7;
    D=-2.27E+20;
    E=4.54E22;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('Dodecane',fuel)
    A =  8.8E-1;
    B = -1.091E3;
    C = -5.03E7;
    D = -5.4874E21;
    E =  1.4959E24;
     
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('Ethanol',fuel)
    A=4.4000E-2;
    B=-5.5700E+1;
    C=-1.2900E+7;
    D=-6.49E+19;
    E=-2.58E+22;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('iso-Octane',fuel)
    A=2.6838E-1;
    B=-3.2388E+2;
    C=-4.6899E+7;
    D=-5.9491E+19;
    E=8.8147E+21;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('THFA',fuel)
    A=1.9766E-1;
    B=-1.8908E+2;
    C=-8.1632E+7;
    D=-5.5857E+20;
    E=8.485E22;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('Heptane',fuel)
    A=2.746E-1;
    B=-2.91E2;
    C=-4.418E+7;
    D=-8.8E19;
    E=1.285E22;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('Decanol',fuel)
    A=4.8836E-1;
    B=-5.8491E+2;
    C=-1.1672E+8;
    D=3.8936E+21;
    E=-1.5303E24;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
   
elseif strcmp('Ethyl-Levulinate',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('Butyl-Levulinate',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('2-MTHF',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
    
elseif strcmp('DMS-T-1.5',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;
elseif strcmp('DMS-T-2',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;

elseif strcmp('DMS-T-5',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;

elseif strcmp('DMS-T-10',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;

elseif strcmp('DNBE',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;

elseif strcmp('Diesel',fuel)
 A=0;
    B=0;
    C=0;
    D=0;
    E=0;
    B2 = A + B/T + C/T^3 + D/T^8 + E/T^9;

end
end




