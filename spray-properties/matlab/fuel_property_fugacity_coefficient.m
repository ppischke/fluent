function [Phi]        = fuel_property_fugacity_coefficient(T,p,fuel);

% [B2] = fuel_property_virial_coefficient(T,fuel);
% B2   = B2/1000; %m^3/mol
% R    = 8.314;   %J/molK
% Phi  = exp(2*B2*p/R/T-log(1+B2*p/R/T));
% %Phi=exp(2*B2*R*T/p*(1+B2*p/R/T)-log(1+B2*p/R/T));

[a,b]   = fuel_property_van_der_vaals_constants(fuel);
R=8.314*1000; %J/kmolK
V=R*T/p; %m3/kmol
Z=V/(V-b)-a/(R*T*V); %Gl. 2.19 Pfennig
Phi = exp(b/(V-b)-log(Z*(1-b/V))-2*a/(R*T*V));

end