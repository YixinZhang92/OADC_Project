clear all; close all; clc;

data = textread('src-2.10');
x = data(:,5);
y = data(:,6);
z = -data(:,7);

outfile = 'CSZ_hypos.txt';
% 
% % plot the input data
% % 
% figure;
% plot3(ys,xs,-zs,'o');
% axis equal;
% grid on;
% title('Input hypocenters');
% xlabel('X km');
% ylabel('Y km');
% zlabel('Z km');


scatter3(x,y,z,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); %(Coordiante (0,0) is at the lower right corner)
grid MINOR; hold on


% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:length(x)
    
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[x(kk) y(kk) z(kk)]);
    
end
