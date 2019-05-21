clear all; close all; clc;

% Method 1
% infile = 'M3.MCU_mincdf_0.80286.txt';
% simul_tag = 'M3.MCUCO.vol0.3';

% % Method 2
infile = 'CSZ_hypos.txt';
simul_tag = 'BigData';

% infile = 'M3.MCUtestBigData_mincdf_0.84274.txt';
% simul_tag = 'BigData';


n_iter =5; uncert = 1.5; % parameters for collapsing algorithm
skim = 'Nichol_etal'; % Choose the skim to use.
             % 'JS' = Jone and Stewart (1997) (Centroid = mean)
             % 'Nichol_etal' = Nicholson et al.(2000) (Centroid = Weighted mean)
                
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

%% Assume that the uncertainty ellipsoid of eac hypocenter is 1km
[x_colap, y_colap, z_colap, dist_to_ref_loc] = calc_collapsed_cloud(xs, ys, zs, uncert, n_iter, skim);

%% Plot original and collapsed hypocenters with the probability functions
for k = 1:n_iter
    
    %if k==n_iter
        % Wrte collapsed hypocenters to file
        outfile = [simul_tag '_uncert_' num2str(uncert) 'km_niter_' num2str(k) '.txt'];
        fid=fopen(outfile,'w');
        for kk=1:length(x_colap(k,:))

            fprintf(fid,'%12.5f %12.5f %12.5f\n',[x_colap(k,kk) y_colap(k,kk) z_colap(k,kk)]);

        end
    %end
    
    % Extract distances moved by each hypocenter at the current iteration
    % relative to its original location. Remove those that have not moved.
    a = dist_to_ref_loc(k,:);
    a(a<0.0001) = [];
    
    FigH = figure('Position', get(0, 'Screensize'));
    %figure
    
%     ax1 = subplot (2,2,1);
    pos1 = [0.05 0.55 0.4 0.4];
    ax1 = subplot('Position',pos1);
    
    scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
    xlabel('Longitude (km)','FontSize',14)
    ylabel('Latitude (km)','FontSize',14)
    zlabel('depth (km)','FontSize',14)
    title(['Earthquake Distribution of CSZ (' simul_tag,')'],'FontSize',18); 
    set(gca, 'fontsize', 14);
    view([-75.9 68.791015625]);
    
    %hold on;
%     ax2 = subplot (2,2,2);
    pos1 = [0.53 0.55 0.4 0.4];
    ax2 = subplot('Position',pos1);
    
    scatter3(x_colap(k,:),y_colap(k,:),z_colap(k,:),'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 .75 .75])
    xlabel('Longitude (km)','FontSize',14)
    ylabel('Latitude (km)','FontSize',14)
    zlabel('depth (km)','FontSize',14)
    %title(['Collapsed Earthquake Distribution ( no of iteration = ', num2str(k),' )'],'FontSize',18); 
    title(['Collapsed Earthquake Distribution (n\_iter = ', num2str(k),')'],'FontSize',18); 
    set(gca, 'fontsize', 14);
    view([-75.9 68.791015625]);
    
    Link = linkprop([ax1, ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    setappdata(gcf, 'StoreTheLink', Link);

    %subplot(2,2,3)
    len = 0.27;
    pos1 = [0.05 0.1 len len];
    subplot('Position',pos1)
    pd1 = createFit_pdf(a);
    title('Probability Density','FontSize',18); 
    set(gca, 'fontsize', 14);
    
    %subplot(2,2,4)
    pos2 = [0.365 0.1 len len];
    subplot('Position',pos2)
    pd2 = createFit_cdf(a);
    title('Cummulative Probability Density','FontSize',18); 
    set(gca, 'fontsize', 14);
    
%     subplot(1,3,3)
    pos2 = [0.7 0.1 len len];
    subplot('Position',pos2)
    pd3 = createFit_pf(a);
    title('Probability Plot','FontSize',18); 
    set(gca, 'fontsize', 14);
    
    % Printing figure to file
    fig_filename = [simul_tag '_uncert_' num2str(uncert) 'km_niter_' num2str(k) '.png'];
    F    = getframe(FigH);
    imwrite(F.cdata, fig_filename, 'png')
    
    a=[];
end
