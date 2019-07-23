clear all; close all; clc;

cdf_max = 0.6;

[xs,ys,zs] = creating_textfile_for_CSZ('CSZ_hypos.txt');
[xs_sort, ys_sort, zs_sort, vol_sort] = calc_tetvol(xs,ys,zs);

figure
pd = makedist('Normal'); p = cdf(pd,vol_sort);
plot(vol_sort,p,'b-','LineWidth',3); grid MINOR; 
xlabel('Volume (km3)','FontSize',18); xlim([0 5]);
ylabel('Cumulative distribution of tetrahedra volumes','FontSize',18)
title('Cumulative Distribution','FontSize',18); 
set(gca, 'fontsize', 18);

xs_declu = xs_sort(p<=cdf_max);
ys_declu = ys_sort(p<=cdf_max);
zs_declu = zs_sort(p<=cdf_max);

xs_unclu = xs_sort(p>cdf_max);
ys_unclu = ys_sort(p>cdf_max);
zs_unclu = zs_sort(p>cdf_max);

%% write synthetic data to outfile
outfile = ['declustered_hypo_' num2str(cdf_max) '.txt'];
fid=fopen(outfile,'w');
for kk=1:length(xs_declu)
    
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_declu(kk) ys_declu(kk) zs_declu(kk)]);
    
end

%% Figures
figure
ax1 = subplot (2,2,1);
scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); 
%set(gca, 'fontsize', 18);
%(Coordiante (0,0) is at the lower right corner)
grid MINOR; hold on

ax2 = subplot (2,2,3);
scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18); 
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title(['Clustered Earthquakes (cdf max = ' num2str(cdf_max) ')'],'FontSize',18); 
grid MINOR; hold on

ax3 = subplot (2,2,2);
scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]); hold on;
scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor','r'); hold on;   
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Clustered and Unclustered Earthquakes','FontSize',18); 
grid MINOR; hold on

ax4 = subplot (2,2,4);
scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor','r')
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Unclustered Earthquakes','FontSize',18); 
grid MINOR; hold on

Link = linkprop([ax1, ax2, ax3, ax4], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);

