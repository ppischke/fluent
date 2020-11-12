function[Heat_Vap]=fuel_property_enthalpy_vaporization(T,fuel)
if strcmp('Butanol',fuel)
    A=6.739E+07;
    B=1.73E-01;
    C=2.915E-01;
    D=0;
    E=0;
    T_c=562; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 74.12;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg
elseif strcmp('Dodecane',fuel)
    A=7.7337*10^7;
    B=4.0681*10^-1;
    C=0;
    D=0;
    E=0;
    T_c=658.00; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 170.34;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg
elseif strcmp('Ethanol',fuel)
    A=5.69E+07;
    B=3.359*10^-1;
    C=-1.025*10^-7;
    D=0;
    E=0;
    T_c=514.00; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 46.07;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg
elseif strcmp('iso-Octane',fuel)
    A=4.7711*10^7;
    B=3.7949*10^-1;
    C=0;
    D=0;
    E=0;
    T_c=543.80; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 114.23;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg    
elseif strcmp('THFA',fuel)
    A=6.4109*10^7;
    B=2.8538*10^-1;
    C=0;
    D=0;
    E=0;
    T_c=639.00; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 102.13;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg   
elseif strcmp('Heptane',fuel)
    A=50014000;
    B=0.38795;
    C=0;
    D=0;
    E=0;
    T_c=540.20; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 100.21;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg 
elseif strcmp('Decanol',fuel)
    A=117500000;
    B=0.65112;
    C=0;
    D=0;
    E=0;
    T_c=687.30; %Kritische Temperatur in K
    T_r=T/T_c;
    M_fuel= 158.28;
    Heat_Vap_mol = A*(1-T_r)^(B + C*T_r + D*T_r^2 + E*T_r^3)/1000; %J/mol
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg 
elseif strcmp('2-MTHF',fuel)
    %Heat_Vap     = -0.000055923*T^4+0.072491*T^3-34.882*T^2+6730*T; %J/kg 
    Heat_Vap     = -5.5923E-5*T^4+0.072491*T^3-34.882*T^2+6730*T; %J/kg  %Kai Leonhardt
elseif strcmp('DNBE',fuel)
    A=59616000;
    B=0.38833;
    T_c=584.1;
    T_r=T/T_c;
    M_fuel= 130.23;
    Heat_Vap_mol = (A*(1-T_r)^B)/1000;
    Heat_Vap     = Heat_Vap_mol/M_fuel*1000; %J/kg 
end
end
