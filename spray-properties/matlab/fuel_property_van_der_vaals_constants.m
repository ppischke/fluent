function [a,b]        = fuel_property_van_der_vaals_constants(fuel);
R=8.314*1000;%J/kmolK
if strcmp('Butanol',fuel)
   a= 20.94/10*1000^2; %Pa m6/kmol2 
   b= 0.1326/1000*1000; %m3/kmol Van der Vaals Volume

elseif strcmp('Dodecane',fuel)
   a= 69.38/10*1000^2; %Pa m6/kmol2 
   b= 0.3758/1000*1000; %m3/kmol Van der Vaals Volume
    
elseif strcmp('Ethanol',fuel)
   a= 12.56/10*1000^2; %Pa m6/kmol2 
   b= 0.0871/1000*1000; %m3/kmol Van der Vaals Volume
    
elseif strcmp('iso-Octane',fuel)
    [T_c]=fuel_property_T_crit(fuel);
    [p_c]=fuel_property_p_crit(fuel);
  
    a=27/64*(R*T_c)^2/p_c;
    b=0.125*R*T_c/p_c;

elseif strcmp('THFA',fuel)
    [T_c]=fuel_property_T_crit(fuel);
    [p_c]=fuel_property_p_crit(fuel);
  
    a=27/64*(R*T_c)^2/p_c;
    b=0.125*R*T_c/p_c;
    
elseif strcmp('Heptane',fuel)
   a= 31.06/10*1000^2; %Pa m6/kmol2 
   b= 0.2049/1000*1000; %m3/kmol Van der Vaals Volume
       
elseif strcmp('Decanol',fuel)
   a= 59.51/10*1000^2; %Pa m6/kmol2 
   b= 0.3086/1000*1000; %m3/kmol Van der Vaals Volume
 
elseif strcmp('2-MTHF',fuel)
    [T_c]=fuel_property_T_crit(fuel);
    [p_c]=fuel_property_p_crit(fuel);
  
    a=27/64*(R*T_c)^2/p_c;
    b=0.125*R*T_c/p_c;
end
end




