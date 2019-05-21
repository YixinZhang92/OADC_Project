clear all; close all; clc;

infile = 'CSZ_hypos.txt';
sigma_x = 0.15;
sigma_y = sigma_x;
sigma_z = 0.35;

%% OPEN HYPOCENTER FILE AND PLOT
fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

org_cat = [xs' ys' zs'];

%% Covariance materix
cov_mat1 = [sigma_x sigma_y sigma_z 0 0 0];
cov_mat = repmat(cov_mat1, N,1);
cov_mat(1,1) = cov_mat(1,1)- 0.001;

%% Applying condensation algorithm
con_cat=catCondens(org_cat,cov_mat,1);

















% % Assume that the uncertainty ellipsoid of eac hypocenter is 1km
% [x_colap, y_colap, z_colap] = calc_collapsed_cloud(xs, ys, zs, 1., 10);
% 
% figure
% 
% ax1 = subplot (2,2,1);
% scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Clustered Earthquake Distribution (from above CDF method)','FontSize',18); 
% 
% hold on;
% ax2 = subplot (2,2,2);
% scatter3(x_colap,y_colap,z_colap,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 1 .0])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Collapsed Earthquake Distribution','FontSize',18); 
% 
% Link = linkprop([ax1, ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
% setappdata(gcf, 'StoreTheLink', Link);
% 
% 
% 
% 
% %% write synthetic data to outfile
% outfile = ['declustered_collapsed_hypo_' num2str(cdf_max) '.txt'];
% fid=fopen(outfile,'w');
% for kk=1:length(xs_declu)
%     
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_declu(kk) ys_declu(kk) zs_declu(kk)]);
%     
% end
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% [xs_sort, ys_sort, zs_sort, vol_sort] = calc_tetvol(x_colap,y_colap,z_colap);
% 
% figure
% pd = makedist('Normal'); p = cdf(pd,vol_sort);
% plot(vol_sort,p,'b-','LineWidth',3); grid MINOR; 
% xlabel('Volume (km3)','FontSize',18)
% ylabel('Cumulative distribution of tetrahedra volumes','FontSize',18)
% title('Cumulative Distribution of collapsed earthquake','FontSize',18); 
% set(gca, 'fontsize', 18);
% 
% xs_declu = xs_sort(p<=cdf_max);
% ys_declu = ys_sort(p<=cdf_max);
% zs_declu = zs_sort(p<=cdf_max);
% 
% xs_unclu = xs_sort(p>cdf_max);
% ys_unclu = ys_sort(p>cdf_max);
% zs_unclu = zs_sort(p>cdf_max);
% 
% %% write synthetic data to outfile
% outfile = ['declustered_collapsed_hypo_' num2str(cdf_max) '.txt'];
% fid=fopen(outfile,'w');
% for kk=1:length(xs_declu)
%     
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_declu(kk) ys_declu(kk) zs_declu(kk)]);
%     
% end
% 
% 
% outfile = ['collapsed_hypo_' num2str(cdf_max) '.txt'];
% fid=fopen(outfile,'w');
% for kk=1:length(x_colap)
%     
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',[x_colap(kk) y_colap(kk) z_colap(kk)]);
%     
% end
% 
% 
% % xs_unclu = xs_sort(p>cdf_max);
% % ys_unclu = ys_sort(p>cdf_max);
% % zs_unclu = zs_sort(p>cdf_max);
% 
% %%
% figure
% 
% ax1 = subplot (2,2,1);
% scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Clustered Earthquake Distribution (from above CDF method)','FontSize',18); 
% 
% hold on;
% ax2 = subplot (2,2,2);
% scatter3(x_colap,y_colap,z_colap,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Collapsed Earthquake Distribution','FontSize',18); 
% 
% ax3 = subplot (2,2,3);
% scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]); hold on;
% scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor','r'); hold on;   
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Clustered and Unclustered Earthquakes (collapsed)','FontSize',18); 
% grid MINOR; hold on
% 
% hold on;
% ax4 = subplot (2,2,4);
% scatter3(xs_declu,ys_declu,zs_declu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 0 1])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title(['Final Clustered Earthquakes (cdf max = ' num2str(cdf_max) ')'],'FontSize',18); 
% 
% 
% Link = linkprop([ax1, ax2, ax3, ax4], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
% setappdata(gcf, 'StoreTheLink', Link);
% 
