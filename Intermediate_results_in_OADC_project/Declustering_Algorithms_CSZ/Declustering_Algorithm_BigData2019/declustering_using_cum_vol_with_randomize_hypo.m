clear all; close all; clc;

cdf_max = 0.85;

[xs,ys,zs] = creating_textfile_for_CSZ('CSZ_hypos.txt');
[xs_rand, ys_rand, zs_rand]=randomize_hypocenters(xs,ys,zs,100);




[xs_sort, ys_sort, zs_sort, vol_sort] = calc_tetvol(xs,ys,zs);
[xsrand_sort, ysrand_sort, zsrand_sort, volrand_sort] = calc_tetvol(xs_rand,ys_rand,zs_rand);



figure
pd = makedist('Normal'); p = cdf(pd,vol_sort);
pd_rand = makedist('Normal'); p_rand = cdf(pd_rand,volrand_sort);
plot(vol_sort,p,'r-'); hold on;
plot(volrand_sort,p_rand,'b-'); 
















xs_declu = xs_sort(p<=cdf_max);
ys_declu = ys_sort(p<=cdf_max);
zs_declu = zs_sort(p<=cdf_max);

xs_unclu = xs_sort(p>cdf_max);
ys_unclu = ys_sort(p>cdf_max);
zs_unclu = zs_sort(p>cdf_max);



figure
ax1 = subplot (2,2,1);
scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); 
%(Coordiante (0,0) is at the lower right corner)
grid MINOR; hold on

ax2 = subplot (2,2,2);
scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Declustered Earthquake Distribution','FontSize',18); 
%(Coordiante (0,0) is at the lower right corner)
grid MINOR; hold on

ax3 = subplot (2,2,3);
scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Unclustered Earthquake Distribution','FontSize',18); 
%(Coordiante (0,0) is at the lower right corner)
grid MINOR; hold on


Link = linkprop([ax1, ax2, ax3], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);

