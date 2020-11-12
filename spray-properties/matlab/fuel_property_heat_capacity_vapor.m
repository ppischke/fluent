function [CPV]= fuel_property_heat_capacity_vapor(T,fuel)

%if fuel==char('Butanol')
if strcmp('Butanol',fuel)
    A=7.4540E4;
    B=2.5907E5;
    C=1.6073E3;
    D=1.7320E5;
    E=7.1240E2;
    M_fuel=74.12;
    CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
elseif strcmp('Dodecane',fuel)
    A=2.1295E5;
    B=6.633E5;
    C=1.7155E3;
    D=4.5161E5;
    E=7.775E2;
    
    M_fuel=170.34;
   CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
elseif strcmp('Ethanol',fuel)
    A=4.9200E4;
    B=1.4577E5;
    C=1.6628E3;
    D=9.3900E4;
    E=7.447E2;
    M_fuel=46.07;
 CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
 
elseif strcmp('iso-Octane',fuel)
    A=1.139E5;
    B=5.286E5;
    C=1.594E3;
    D=3.351E5;
    E=6.7794E2;
    M_fuel=114.23;
   CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
      
elseif strcmp('THFA',fuel)
    A=2.6823E5;
    B=1.1670E5;
    C=1.686E3;
    D=-3.52E5;
    E=2.452E2;
    M_fuel=102.13;
   CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
  
elseif strcmp('Heptane',fuel)
    A=1.2015E5;
    B=4.001E5;
    C=1.6766E3;
    D=2.74E5;
    E=7.564E2;
    M_fuel=100.20;
   CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
     
elseif strcmp('Decanol',fuel)
    A=1.6984E5;
    B=5.392E5;
    C=1.568E3;
    D=3.938E5;
    E=7.205E2;
    M_fuel=158.28;
   CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
     
elseif strcmp('DNBE',fuel)
    A=1.6122E5;
    B=4.4777E+5;
    C=1.6831E+3;
    D=2.9180E+5;
    E=7.8160E+2;
    M_fuel=130.23;
   CPV = (A+B*((C/T)/sinh(C/T))^2+D*((E/T)/cosh(E/T))^2)/M_fuel;
     
elseif strcmp('2-MTHF',fuel) %Kai Leonhardt Excel File
    A=4.6801E-08;
    B=-2.5317E-04;
    C=5.0096E-01;
    D=-1.8322E+01;
    M_fuel=86.1323;%g/mol
    CPV = (A*T^3+B*T^2+C*T+D)/(M_fuel/1000);
end
end




