number = 24;

rho    = 680;
temp   = 300;
flow   = 0.0351;
vel    = 150;
delta  = 0.00e-3;
rad    = 2.00e-3;
rot    = 22.5;
cone   = 47.5 + 0.0*randn(1,number);
diam   = flow / 2/pi/rad / rho / vel;

phi    = ((1:number)-(1+number)/2)/(number/2) * rot/180*pi;

x_x    = (delta*sin(cone/180*pi) + rad) .* cos(phi);
x_y    = (delta*sin(cone/180*pi) + rad) .* sin(phi);
x_z    = (delta*cos(cone/180*pi)) .* ones(1,number);

u_x    = vel * sin(cone/180*pi) .* cos(phi);
u_y    = vel * sin(cone/180*pi) .* sin(phi);
u_z    = vel * cos(cone/180*pi) .* ones(1,number);

f = fopen(sprintf('u:\\Fluent12\\spray-udf\\trunk\\nozzle_gdi.inj'),'w');

for i=1:number
    
    fprintf(f,'((%e %e %e %e %e %e %e %e %e) injection:%d)\n',...
        x_x(i),x_y(i),x_z(i),...
        u_x(i),u_y(i),u_z(i),diam,temp,flow*rot/180/number,i);
    
end

fclose(f);



f = fopen(sprintf('u:\\Fluent12\\spray-udf\\jobs\\demo-gdi\\vd\\n%.2d\\nozzle_gdi.inj',number/12),'w');

for i=1:number
    
    fprintf(f,'((%e %e %e %e %e %e %e %e %e) injection:%d)\n',...
        x_x(i),x_y(i),x_z(i),...
        u_x(i),u_y(i),u_z(i),diam,temp,flow*rot/180/number,i);
    
end

fclose(f);
