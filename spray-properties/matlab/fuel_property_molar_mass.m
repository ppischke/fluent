function[M_fuel]=fuel_property_molar_mass(fuel)
if strcmp('Butanol',fuel)

    M_fuel=74.1214; %g/mol

elseif strcmp('Ethanol',fuel)

    M_fuel=46.07; %g/mol

elseif strcmp('THFA',fuel)

    M_fuel=102.13; 

elseif strcmp('iso-Octane',fuel)

    M_fuel= 114.23; %g/mol

elseif strcmp('Dodecane',fuel)

    M_fuel=170.34; %g/mol

elseif strcmp('Heptane',fuel)

    M_fuel=100.21; %g/mol

elseif strcmp('Decanol',fuel)

    M_fuel=158.28; %g/mol

elseif strcmp('2-MTHF',fuel)
    
    M_fuel=86.1323;%g/mol

elseif strcmp('DNBE',fuel)

    M_fuel=130.23;

end
end


