clear all
fuel={char('Decanol'), char('Butanol'),char('Ethanol'),char('Dodecane'),char('THFA'),char('Heptane'),char('iso-Octane'),char('DNBE'),char('2-MTHF')};
%fuel={char('Decanol'), char('Butanol'),char('Ethanol')};


    
%% Figuredesign
gcf = figure;
set(gcf, 'InvertHardcopy','on', 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters');
%set(gcf, 'PaperType', 'A4','PaperOrientation','portrait');
set(gcf, 'PaperSize', [11.9 5.5]);
%set(gcf, 'PaperSize', [13.5 8.5]);
set(gcf, 'PaperPosition', [1 -0.2 10 6]);
axes1 = axes('Parent',gcf,'FontSize',8);
%axes1 = axes('Parent',gcf,'Position',[0.2 0.1 0.6 0.5],'FontSize',12);
markersize=3;
linienstaerke = 0.5;
RGB1=[0, 0, 0];
RGB1f=[1 1 1];
RGB1e=[0 0 0];
%xlim([0 1.5]); ylim([0 40]); 
box('on');
xlabel({'Temperatur / K',''},'FontName','Arial','FontSize',8);
ylabel({'Verdampfungsenthalpie / kJ/kg',''},'FontName','Arial','FontSize',8);
%ylabel({'Dampfdruck / MPa',''},'FontName','Arial','FontSize',8);

%%
for i=1:length(fuel)
    
[T_c]=fuel_property_T_crit(fuel{i});
[T_cref]=fuel_property_T_crit(fuel{1});
[pV_c]=fuel_property_vapor_pressure(T_c,fuel{i});
[pV_cref]=fuel_property_vapor_pressure(T_c,fuel{1});
T=[300:10:T_c];
T_ref=[300:10:T_cref];
clear Heat_Vao_T
clear pV_T
    for j=1:length(T_ref)
    pVref=fuel_property_vapor_pressure(T_ref(j),fuel{1});
    pVref_T(j)=pVref/1000/1000;
    end
    
    for j=1:length(T)
    [Heat_Vap]=fuel_property_enthalpy_vaporization(T(j),fuel{i});
    [pV]=fuel_property_vapor_pressure(T(j),fuel{i});
    Heat_Vao_T(j)=Heat_Vap/1000;
    pV_T(j)=pV/1000/1000;
    end
    j=j+1;
    Heat_Vao_T(j)=fuel_property_enthalpy_vaporization(T_c,fuel{i})/1000;
    pV_T(j)=fuel_property_vapor_pressure(T_c,fuel{i})/1000/1000;
    pVref_T(j)=fuel_property_vapor_pressure(T_cref,fuel{1})/1000/1000;
    T(j)=T_c;
    hold on
 %T=T./T_c;
 %pV_T=pV_T./pV_c;
 pVref_T=pVref_T./pV_cref;
 
 %pV_T=pV_T./pVref_T;
 
if strcmp('Butanol',fuel{i})
    exp_plot_Butanol=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','o','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    %exp_plot_Butanol=plot(T,pV_T,'MarkerSize',markersize,'LineWidth',1,'Marker','o','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    %exp_plot_Butanol=plot(T,pV_T,'MarkerSize',markersize,'Marker','o','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_Butanol={[char(fuel{i})]}; 
elseif strcmp('Ethanol',fuel{i})
    exp_plot_Ethanol=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','<','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    %exp_plot_Ethanol=plot(T,pV_T,'MarkerSize',markersize,'Marker','<','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_Ethanol={[char(fuel{i})]}; 
elseif strcmp('THFA',fuel{i})
    exp_plot_THFA=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','>','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    %exp_plot_THFA=plot(T,pV_T,'MarkerSize',markersize,'Marker','>','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_THFA={[char(fuel{i})]}; 
elseif strcmp('2-MTHF',fuel{i})
    exp_plot_2_MTHF=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','^','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
   % exp_plot_2_MTHF=plot(T,pV_T,'MarkerSize',markersize,'Marker','^','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_2_MTHF={[char(fuel{i})]}; 
elseif strcmp('DNBE',fuel{i})
    exp_plot_DNBE=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','v','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
   % exp_plot_DNBE=plot(T,pV_T,'MarkerSize',markersize,'Marker','v','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_DNBE={[char(fuel{i})]}; 
elseif strcmp('iso-Octane',fuel{i})
    exp_plot_iso_Octane=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','s','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
   % exp_plot_iso_Octane=plot(T,pV_T,'MarkerSize',markersize,'Marker','s','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_iso_Octane={[char(fuel{i})]}; 
elseif strcmp('Heptane',fuel{i})
    exp_plot_Heptane=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','h','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    %exp_plot_Heptane=plot(T,pV_T,'MarkerSize',markersize,'Marker','h','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_Heptane={[char(fuel{i})]}; 
elseif strcmp('Dodecane',fuel{i})
    exp_plot_Dodecane=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','x','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
   % exp_plot_Dodecane=plot(T,pV_T,'MarkerSize',markersize,'Marker','x','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_Dodecane={[char(fuel{i})]}; 
elseif strcmp('Decanol',fuel{i})
    exp_plot_Decanol=plot(T,Heat_Vao_T,'MarkerSize',markersize,'Marker','p','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
   % exp_plot_Decanol=plot(T,pV_T,'MarkerSize',markersize,'Marker','p','Color',RGB1e,'MarkerFaceColor',RGB1f,'MarkerEdgeColor',RGB1e);
    legendname_Decanol={[char(fuel{i})]}; 
end
  

end
%legend([exp_plot_Butanol,exp_plot_Ethanol,exp_plot_THFA,exp_plot_2_MTHF,exp_plot_DNBE,exp_plot_iso_Octane,exp_plot_Heptane,exp_plot_Dodecane,exp_plot_Decanol], [legendname_Butanol,legendname_Ethanol,legendname_THFA,legendname_2_MTHF,legendname_DNBE,legendname_iso_Octane,legendname_Heptane,legendname_Dodecane,legendname_Decanol],'Location','NorthWest');
legend([exp_plot_Butanol,exp_plot_Ethanol,exp_plot_THFA,exp_plot_2_MTHF,exp_plot_DNBE,exp_plot_iso_Octane,exp_plot_Heptane,exp_plot_Dodecane,exp_plot_Decanol], [legendname_Butanol,legendname_Ethanol,legendname_THFA,legendname_2_MTHF,legendname_DNBE,legendname_iso_Octane,legendname_Heptane,legendname_Dodecane,legendname_Decanol],'Location','NorthEast');

%% Bildspeicherung
%Dateiname=char(fnames_split(5));
%Dateiname=Dateiname(1:length(Dateiname)-4);
print(gcf, '-dpdf', '-r300', ['Z:\Ergebnisse\figures\Verdampfungsenthalpie.pdf'])
%print(gcf, '-dpdf', '-r300', ['Z:\Ergebnisse\figures\Dampfdruck.pdf'])