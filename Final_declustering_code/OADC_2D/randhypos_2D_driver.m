clear all; close all; clc;

no_simul = 3;

if no_simul==1
    % Simulation 1 - 1 line
    L = 10; err_av=0.2; nhypos= 200; strike= -35; rt = [0 0];
    [rxp, ryp] = rand_hypos_2D(L,err_av,nhypos,strike,rt);

elseif no_simul==2
    % Simulation 2 - 3 lines
    L1 = 10; err_av1=0.2; nhypos1= 200; strike1= -35; rt1 = [0 5];
    [rxp1, ryp1] = rand_hypos_2D(L1,err_av1,nhypos1,strike1,rt1); 
    
    L2 = 10; err_av2=0.2; nhypos2= 200; strike2= -35; rt2 = [3.3 5.5];
    [rxp2, ryp2] = rand_hypos_2D(L2,err_av2,nhypos2,strike2,rt2); 
    
    L3 = 10; err_av3=0.2; nhypos3= 200; strike3= -35; rt3 = [4.5 6];
    [rxp3, ryp3] = rand_hypos_2D(L3,err_av3,nhypos3,strike3,rt3); 

    nhypos = nhypos1+nhypos2+nhypos3;
    rxp = [rxp1 rxp2 rxp3];
    ryp = [ryp1 ryp2 ryp3];

elseif no_simul==3
    % Simulation 3 - 4 cross lines
    L1 = 10; err_av1=0.2; nhypos1= 200; strike1= -35; rt1 = [0 5];
    [rxp1, ryp1] = rand_hypos_2D(L1,err_av1,nhypos1,strike1,rt1); 
    
    L2 = 10; err_av2=0.2; nhypos2= 200; strike2= -35; rt2 = [3 5];
    [rxp2, ryp2] = rand_hypos_2D(L2,err_av2,nhypos2,strike2,rt2); 
    
    L3 = 10; err_av3=0.2; nhypos3= 200; strike3= 35; rt3 = [0 5];
    [rxp3, ryp3] = rand_hypos_2D(L3,err_av3,nhypos3,strike3,rt3); 
    
    L4 = 10; err_av4=0.2; nhypos4= 200; strike4= 35; rt4 = [3 5];
    [rxp4, ryp4] = rand_hypos_2D(L4,err_av4,nhypos4,strike4,rt4); 

    nhypos = nhypos1+nhypos2+nhypos3+nhypos4;
    rxp = [rxp1 rxp2 rxp3 rxp4];
    ryp = [ryp1 ryp2 ryp3 ryp4];
    
elseif no_simul==4
    % Simulation 4 - 7 lines ('real')
    L1 = 10; err_av1=0.2; nhypos1= 200; strike1= -55; rt1 = [0 0];
    [rxp1, ryp1] = rand_hypos_2D(L1,err_av1,nhypos1,strike1,rt1); 
    
    L2 = 7; err_av2=0.2; nhypos2= 200; strike2= -75; rt2 = [1 3.5];
    [rxp2, ryp2] = rand_hypos_2D(L2,err_av2,nhypos2,strike2,rt2); 
    
    L3 = 10; err_av3=0.2; nhypos3= 200; strike3= -75; rt3 = [3 4.5];
    [rxp3, ryp3] = rand_hypos_2D(L3,err_av3,nhypos3,strike3,rt3); 
    
    L4 = 3; err_av4=0.2; nhypos4= 100; strike4= -55; rt4 = [4 -1.5];
    [rxp4, ryp4] = rand_hypos_2D(L4,err_av4,nhypos4,strike4,rt4); 

    L5 = 3; err_av5=0.2; nhypos5= 100; strike5= 45; rt5 = [7 -1];
    [rxp5, ryp5] = rand_hypos_2D(L5,err_av5,nhypos5,strike5,rt5); 
    
    L6 = 2; err_av6=0.2; nhypos6= 100; strike6= -50; rt6 = [6.5 1.5];
    [rxp6, ryp6] = rand_hypos_2D(L6,err_av6,nhypos6,strike6,rt6); 
    
    L7 = 3; err_av7=0.2; nhypos7= 100; strike7= 65; rt7 = [4 1.5];
    [rxp7, ryp7] = rand_hypos_2D(L7,err_av7,nhypos7,strike7,rt7);
    
    nhypos = nhypos1+nhypos2+nhypos3+nhypos4+nhypos5+nhypos6+nhypos7;
    rxp = [rxp1 rxp2 rxp3 rxp4 rxp5 rxp6 rxp7];
    ryp = [ryp1 ryp2 ryp3 ryp4 ryp5 ryp6 ryp7];       
end

%figure;
Fig1 = figure('Name','Projected hypocenters','Position', get(0, 'Screensize'));
plot(rxp,ryp,'o');
axis equal;
title('Final translation of hypocenters');
xlabel('X km');
ylabel('Y km');
grid on;

% write synthetic data to outfile
outfile = ['randhypos.2D.simul.' num2str(no_simul) '.txt'];

fid=fopen(outfile,'w');
for kk=1:nhypos
    
    fprintf(fid,'%12.5f %12.5f\n',[rxp(kk) ryp(kk)]);
    
end
