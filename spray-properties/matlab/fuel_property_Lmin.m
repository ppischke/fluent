function [Lmin]= fuel_property_Lmin(fuel)

if strcmp('Butanol',fuel)
    Lmin= 11.16;
elseif strcmp('Dodecane',fuel)
    Lmin= 14.98 ;
elseif strcmp('Ethanol',fuel)
     Lmin=  8.98;
elseif strcmp('iso-Octane',fuel)
    Lmin=  15.09;
elseif strcmp('THFA',fuel)
     Lmin=  8.78;
elseif strcmp('Heptane',fuel)
     Lmin=  15.14;
elseif strcmp('Decanol',fuel)
   Lmin=  13.07;  
elseif strcmp('Ethyl-Levulinate',fuel)
  Lmin= 8.13; 
elseif strcmp('Butyl-Levulinate',fuel)
  Lmin=  9.21; 
elseif strcmp('2-MTHF',fuel)
  Lmin=  11.21; 
elseif strcmp('DMS-T-1.5',fuel)
  Lmin=NaN;
elseif strcmp('DMS-T-2',fuel)
  Lmin= NaN;
elseif strcmp('DMS-T-5',fuel)
  Lmin= NaN;
elseif strcmp('DMS-T-10',fuel)
  Lmin= NaN;
elseif strcmp('DNBE',fuel)
   Lmin= 12.71;
elseif strcmp('Diesel',fuel)
   Lmin= 14.5;  
end
end




