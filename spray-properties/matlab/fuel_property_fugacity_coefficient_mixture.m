function [Phi_mix]        = fuel_property_fugacity_coefficient_mixture(T,p,fuel,yi);

[B2]        = fuel_property_virial_coefficient(T,fuel);
[a_i,b_i]   = fuel_property_van_der_vaals_constants(fuel);

a_O = 1.378/10*1000^2; %Pa m6/kmol2 Sauerstoff O2
b_O = 0.03183/1000*1000; %m3/kmol Van der Vaals Volume Sauerstoff O2

a_N = 1.408/10*1000^2; %Pa m6/kmol2 Stickstoff N2
b_N = 0.03913/1000*1000; %m3/kmol Van der Vaals Volume Stickstoff N2

yi_i = yi;
yi_O = (1-yi)*0.2095;
yi_N = (1-yi)*0.7808;

a = yi_i*yi_i*a_i + yi_i*yi_O*(a_i*a_O)^(0.5) + yi_i*yi_N*(a_i*a_N)^(0.5) + yi_O*yi_i*(a_O*a_i)^(0.5) + yi_O*yi_O*a_O + yi_O*yi_N*(a_O*a_N)^(0.5) + yi_N*yi_i*(a_N*a_i)^(0.5) + yi_N*yi_O*(a_N*a_O)^(0.5) + yi_N*yi_N*a_N;
b = yi * b_i + yi_N * b_N + yi_O * b_O;

R=8.314*1000; %J/kmolK
V=R*T/p; %m3/kmol

Z=V/(V-b)-a/(R*T*V); %Gl. 2.19 Pfennig

Phi_mix = exp(log(V/(V-b))+b_i/(V-b)-2/R/T/V*(yi_i * a_i+ yi_O *(a_i * a_O)^0.5 + yi_N * (a_i * a_N)^0.5)-log(Z));  %Gl. 4.64 Pfennig

end