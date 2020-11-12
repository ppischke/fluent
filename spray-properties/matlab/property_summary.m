%fuel={char('THFA'),char('Ethyl-Levulinate'),char('Butyl-Levulinate'),char('DMS-T10'),char('DMS-T5'),char('DMS-T2'),char('DMS-T1.5'),char('2-MTHF'),char('Diesel'),char('Decanol'), char('Butanol'),char('Ethanol'),char('DNBE'),char('Dodecane'),char('iso-Octane'),char('Heptane')};
fuel={char('Butanol')};
%i=16;
%fuel{i}
T=300;

[M_fuel]=fuel_property_molar_mass(fuel)
[density]        = fuel_property_density(T,fuel)
[viscosity]= fuel_property_viscosity(T,fuel)*1000
[surface_tension]= fuel_property_surface_tension(T,fuel)*1000


%[M_fuel]=fuel_property_molar_mass(fuel{i})
%[density]        = fuel_property_density(T,fuel{i})
%[viscosity]= fuel_property_viscosity(T,fuel{i})*1000
%[surface_tension]= fuel_property_surface_tension(T,fuel{i})*1000



%[T_c]=fuel_property_T_crit(fuel{i})
%[Heat_Vap]=fuel_property_enthalpy_vaporization(T,fuel{i})
%[pV]=fuel_property_vapor_pressure(T,fuel{i})
%[CP]= fuel_property_heat_capacity(T,fuel{i})
D0=109*10^-6;%m
L0=10*D0;
pU=1*10^5;%Pa
pinje = [20 30 40 50 60 70]
for p=1:6
    pinj=pinje(p)*10^6;
[v0] = outlet_velocity(D0,L0,pinj,pU,density,viscosity/1000)

Re(p) = D0*density*v0/(viscosity/1000);
We(p) = D0*density*v0^2/(surface_tension/1000);
end
We
