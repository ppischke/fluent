function[T_c]=fuel_property_T_crit(fuel)
if strcmp('Butanol',fuel)
    
    T_c=562; %Kritische Temperatur in K
  
elseif strcmp('Dodecane',fuel)
  
    T_c=658.00; %Kritische Temperatur in K
  
elseif strcmp('Ethanol',fuel)
  
    T_c=514.00; %Kritische Temperatur in K
 
elseif strcmp('iso-Octane',fuel)

    T_c=543.80; %Kritische Temperatur in K
 
elseif strcmp('THFA',fuel)
  
    T_c=639.00; %Kritische Temperatur in K
  
elseif strcmp('Heptane',fuel)
 
    T_c=540.20; %Kritische Temperatur in K
  
elseif strcmp('Decanol',fuel)
  
    T_c=687.30; %Kritische Temperatur in K
   
elseif strcmp('2-MTHF',fuel)
  
    T_c=537;
    
elseif strcmp('DNBE',fuel)
   
    T_c=584.1;

end
end
