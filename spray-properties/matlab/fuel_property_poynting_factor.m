function [F]        = fuel_property_poynting_factor(T,pV,pU,fuel);

[M_fuel]  = fuel_property_molar_mass(fuel); %kg/kmol
[density] = fuel_property_density(T,fuel);%kg/m^3

V = M_fuel/density/1000; %m^3/mol
R=8.314; %J/molK

F=exp(V/R/T*(pU-pV));
end