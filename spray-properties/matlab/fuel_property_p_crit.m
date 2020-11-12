function[p_c]=fuel_property_p_crit(fuel)
if strcmp('Butanol',fuel)
    
    %p_c=562; %Kritische Temperatur in K
  
elseif strcmp('Dodecane',fuel)
  
    %p_c=658.00; %Kritische Temperatur in K
  
elseif strcmp('Ethanol',fuel)
  
    p_c=6.148E6;%Pa
 
elseif strcmp('iso-Octane',fuel)

    p_c=2.5700E6;%Pa
 
elseif strcmp('THFA',fuel)
  
    p_c=4.66E6; %Pa
  
elseif strcmp('Heptane',fuel)
 
%    p_c=540.20; %Kritische Temperatur in K
  
elseif strcmp('Decanol',fuel)
  
 %   p_c=687.30; %Kritische Temperatur in K
   
elseif strcmp('2-MTHF',fuel)
  
   p_c=37.58*10^5;
    
elseif strcmp('DNBE',fuel)
   
   % p_c=584.1;

end
end
