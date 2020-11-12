function[pV]=fuel_property_vapor_pressure(T,fuel)
if strcmp('Butanol',fuel)
    A=93.173;
    B=-9185.90;
    C=-9.7464;
    D=4.7796E-18;
    E=6.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('Dodecane',fuel)
    A=137.47;
    B=-11976;
    C=-16.698;
    D=0.0000080906;
    E=2.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('Ethanol',fuel)
    A=74.475;
    B=-7164.3;
    C=-7.327;
    D=0.000003134;
    E=2.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('iso-Octane',fuel)
    A=84.912;
    B=-6722.20;
    C=-9.5157;
    D=0.0000072244;
    E=2.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('THFA',fuel)
    A=150.59;
    B=-11574.00;
    C=-19.03;
    D=0.000014141;
    E=2.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('Heptane',fuel)
    A=87.829;
    B=-6996.40;
    C=-9.8802;
    D=0.0000072099;
    E=2.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('Decanol',fuel)
    A=250.59;
    B=-19169.00;
    C=-32.90;
    D=0.000014627;
    E=2.00;
    pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
elseif strcmp('2-MTHF',fuel)
    %pV = 3.6837E-7*T^4-3.2515E-4*T^3+0.10663*T^2-15.332*T+812.18; %Dampfdruck in Pa %Kai Leonhardt
    pV = (3.6837E-10*T^4-3.2515E-7*T^3+1.0663E-4*T^2-1.5332E-2*T+8.1218E-1)*10^6; %Dampfdruck in Pa %Kai Leonhardt
elseif strcmp('DNBE',fuel)
   A=72.227;
   B=-7537.6;
   C=-7.0596;
   D=9.1442E-18;
   E=6;
   pV = exp(A+B/T+C*log(T)+D*T^E); %Dampfdruck in Pa
end
end


