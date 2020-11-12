function [density_air]=air_density(TU,pU)
RL=287.058;%spezifische Gaskonstante

density_air=pU/(RL*TU);%kg/m^3

%density_air=25;%kg/m^3
%density_air=0.0001*TU^2-0.2073*TU+106.78;%kg/m^3
end