number = 4;

vel    = 380;
temp   = 300;
flow   = 0.0045;
diam   = 140e-6;

x_x    = zeros(1,number);
x_y    = zeros(1,number);
x_z    = zeros(1,number);

u_x    = zeros(1,number);
u_y    = zeros(1,number);
u_z    = ones(1,number)*vel;

f = fopen(sprintf('nozzle_diesel.inj'),'w');

for i=1:number
    
    fprintf(f,'((%e %e %e %e %e %e %e %e %e) injection:%d)\n',...
        x_x(i),x_y(i),x_z(i),...
        u_x(i),u_y(i),u_z(i),diam,temp,flow/number,i);
    
end

fclose(f);



f = fopen(sprintf('nozzle_diesel.inj'),'w');

for i=1:number
    
    fprintf(f,'((%e %e %e %e %e %e %e %e %e) injection:%d)\n',...
        x_x(i),x_y(i),x_z(i),...
        u_x(i),u_y(i),u_z(i),diam,temp,flow/number,i);
    
end

fclose(f);
