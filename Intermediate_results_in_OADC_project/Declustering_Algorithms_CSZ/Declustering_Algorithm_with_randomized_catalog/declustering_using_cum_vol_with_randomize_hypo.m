clear all; close all; clc;

% infile = 'CSZ_hypos.txt';
% simul_tag = 'M3.MCUtestBigData';

infile = 'BigData_uncert_1.5km_niter_5.txt';
simul_tag = 'Collapsed.first';

PROB = 0.05;
factor_of_events = 1.0;
print_clus_t0_file =1;

%% Import natural catalog
%[xs,ys,zs] = creating_textfile_for_CSZ('CSZ_hypos.txt');
fid=fopen(infile,'r');
[data,count]=fscanf(fid,'%g %g %g',[3,inf]);
fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

N = length(xs);

figure;
h1=subplot (1,2,1);
scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18); ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
xlim([min(xs)-50 max(xs)+50]);
ylim([min(ys)-50 max(ys)+50]);
title(['Earthquake Distribution (' simul_tag, ')'],'FontSize',18); 
grid MINOR; hold on; view(0,90) % view vertically downward to remve depth effect
set(gca, 'fontsize', 18);

%%  Pick the polygon
[x_poly,y_poly]=getline(h1,'closed'); hold on
z_poly = zeros(length(x_poly),1);

plot(x_poly,y_poly,'r-');

%% Generating ramdomized catalog
h2=subplot (1,2,2);
[xr,yr,zr]=randomcat(xs(1:round(factor_of_events*N)),ys(1:round(factor_of_events*N)),...
    zs(1:round(factor_of_events*N)),x_poly,y_poly,z_poly);

scatter3(xr,yr,zr,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18); ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Randomized Catalog','FontSize',18); 
grid MINOR; hold on; view(0,90) % view vertically downward to remve depth effect
set(gca, 'fontsize', 18);
hold on;
plot(x_poly,y_poly,'r-');

%% Determine tetrahedra volumes for both natual and randomized catalogs
[xs_sort, ys_sort, zs_sort, vol_sort] = calc_tetvol(xs,ys,zs);
[xsrand_sort, ysrand_sort, zsrand_sort, volrand_sort] = calc_tetvol(xr',yr',zr');

% vol_sort(vol_sort<0.0001) = 0.0001;
% volrand_sort(volrand_sort<0.0001) = 0.0001;

%**************************************************************************
%  Compute the normalized cumulative probability density function of the
%  random catalog tetrahedra volumes
[cdf_rnd,vr]=ecdf(volrand_sort);

nkh=length(vr);
% for kh=1:nkh;
%     fprintf('cdf_rnd %g  vr %g\n',cdf_rnd(kh),vr(kh));
% end

%**************************************************************************
%  Compute the normalized cumulative probability density function of the
%  original catalog tetrahedra volumes0
[cdf_cat,vs]=ecdf(vol_sort);
 
%**************************************************************************
%  Calculate the volume of tetrahedra at the 5% probability level (0.05) 
%  from the random catalog
for kh=1:nkh-1
    if cdf_rnd(kh) < PROB && cdf_rnd(kh+1) >= PROB
        ncalc=kh;
    end
end
ncalc;
V05=vr(ncalc);
fprintf('ncalc = %i, Volume at 5 percent probability (V05, km cubed) = %g\n',ncalc,V05);

%V05=0.4;
%**************************************************************************
%  Calculate N(V(0.05)) from original catalog
nks=length(vs);
for ks=1:nks-1
    if vs(ks) < V05 && vs(ks+1) >= V05
        ncat=ks;
    end
end
ncat;
NV05_cat=cdf_cat(ncat);
fprintf('ncat = %i, Probability at target volume (NV05) = %g\n',ncalc,NV05_cat);

% plot the cdfs with determined values
figure;
%hold on;
%plot(log(vr),cdf_rnd,'-k',log(vs),cdf_cat);
semilogx(vr,cdf_rnd,'-k','LineWidth',2);hold on;
semilogx(vs,cdf_cat,'-b','LineWidth',2);
%plot([log(V05) log(V05)],[0 1],'-r');
semilogx([V05 V05],[0 1],'-r','LineWidth',2);
hold off;

title('Cumulative Distribution Functions');
xlabel('volume (km^3)');
ylabel('Probability');
legend('Random','Catalog');
grid MINOR;
set(gca, 'fontsize', 18);

%%
xs_declu = xs_sort(vol_sort <= V05);
ys_declu = ys_sort(vol_sort <= V05);
zs_declu = zs_sort(vol_sort <= V05);

xs_unclu = xs_sort(vol_sort > V05);
ys_unclu = ys_sort(vol_sort > V05);
zs_unclu = zs_sort(vol_sort > V05);

%% write synthetic data to outfile
if print_clus_t0_file == 1
    outfile = [simul_tag '_mincdf_' num2str(NV05_cat) '.txt'];
    fid=fopen(outfile,'w');
    for kk=1:length(xs_declu)

        fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_declu(kk) ys_declu(kk) zs_declu(kk)]);

    end
end


%% Figures
figure

pos1 = [0.05 0.1 0.4 0.8];
ax1 = subplot('Position',pos1);

scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title(['Earthquake Distribution (' simul_tag, ')'],'FontSize',18); 
grid MINOR; hold on
view([-75.9 68.791015625]);


pos1 = [0.53 0.1 0.4 0.8];
ax2 = subplot('Position',pos1);

scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18); 
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title(['Clustered Earthquakes (min cdf = ' num2str(NV05_cat),')'],'FontSize',18); 
grid MINOR;
view([-75.9 68.791015625]);
    
   Link = linkprop([ax1, ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);

%% Figures
figure

pos1 = [0.05 0.1 0.4 0.8];
ax3 = subplot('Position',pos1);

scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]); hold on;
scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor','r'); hold on;   
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Clustered and Unclustered Earthquakes','FontSize',18); 
grid MINOR; hold on
view([-75.9 68.791015625]);


pos1 = [0.53 0.1 0.4 0.8];
ax4 = subplot('Position',pos1);

scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor','r')
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Unclustered Earthquakes','FontSize',18); 
grid MINOR; hold on
view([-75.9 68.791015625]);
    
Link = linkprop([ax3, ax4], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);


%% Figures
% figure
% ax1 = subplot (2,2,1);
% scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title(['Earthquake Distribution (' simul_tag, ')'],'FontSize',18); 
% grid MINOR; hold on
% 
% ax2 = subplot (2,2,3);
% scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18); 
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title(['Clustered Earthquakes (min cdf = ' num2str(NV05_cat),')'],'FontSize',18); 
% grid MINOR; hold on
% 
% ax3 = subplot (2,2,2);
% scatter3(xs_declu, ys_declu, zs_declu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]); hold on;
% scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor','r'); hold on;   
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Clustered and Unclustered Earthquakes','FontSize',18); 
% grid MINOR; hold on
% 
% ax4 = subplot (2,2,4);
% scatter3(xs_unclu, ys_unclu, zs_unclu,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor','r')
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Unclustered Earthquakes','FontSize',18); 
% grid MINOR; hold on
% 
% Link = linkprop([ax1, ax2, ax3, ax4], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
% setappdata(gcf, 'StoreTheLink', Link);
