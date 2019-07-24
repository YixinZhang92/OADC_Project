function [x_colap, y_colap, z_colap, nSD_moved] = decluster_w_collapsing_JS_EG_3D...
    (infile, simul_tag_col, n_iter, nSD, sigma_x, sigma_y, sigma_z, skim)

%clear all; close all; clc;

% Method 1
% infile = 'M3.MCU_mincdf_0.80286.txt';
% simul_tag = 'M3.MCUCO.vol0.3';

% % % Method 2
% infile = 'CSZ_hypos.txt';
% simul_tag = 'BigData';
% 
% % infile = 'M3.MCUtestBigData_mincdf_0.84274.txt';
% % simul_tag = 'BigData';
% 
% 
% n_iter =1; uncert = 1.5; % parameters for collapsing algorithm
% skim = 'Nichol_etal'; % Choose the skim to use.
%              % 'JS' = Jone and Stewart (1997) (Centroid = mean)
%              % 'Nichol_etal' = Nicholson et al.(2000) (Centroid = Weighted mean)
%            
             
%% OPEN HYPOCENTER FILE AND PLOT
fid=fopen(infile,'r');
[data,~]=fscanf(fid,'%g %g %g',[3,inf]);
fclose(fid);

data=data';
[N,~]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);


%% Collapsing the hypocenters
[x_colap, y_colap, z_colap, nSD_moved] = calc_collapsed_cloud_3D...
    (xs, ys, zs, nSD, sigma_x, sigma_y, sigma_z, n_iter, skim);


%% Print collapsed hypocenters to file
for k = 1:n_iter
    
    %if k==n_iter
        % Wrte collapsed hypocenters to file
        outfile = [simul_tag_col '_niter_' num2str(k) '.txt'];%'_uncert_' num2str(uncert)
        fid=fopen(outfile,'w');
        for kk=1:length(x_colap(k,:))

            fprintf(fid,'%12.5f %12.5f %12.5f\n',[x_colap(k,kk) y_colap(k,kk) z_colap(k,kk)]);

        end
        
        fclose(fid);

end


%% Create theoretical chi distribution
[x,chiscaled] = create_theor_chi_dist(N,3);

% stdn = sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2);
% chiscaled = chiscaled*sqrt(N)/stdn;
chiscaled = chiscaled*100;

%% figures
FigH = figure('Position', get(0, 'Screensize'));

subplot(1,2,1)
scatter3(data(:,1), data(:,2), data(:,3),'filled','MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75]); grid MINOR;
%xlim([-4 4]); ylim([-4 4]); zlim([-4 4]);
xlabel('Longitude (km)','FontSize',14)
ylabel('Latitude (km)','FontSize',14)
zlabel('depth (km)','FontSize',14)
title(['Earthquake Distribution of CSZ (' simul_tag_col,')'],'FontSize',18); 
set(gca, 'fontsize', 14); view([-75.9 68.791015625]);
    
subplot(1,2,2)
plot(x,chiscaled); shg
ylim([0 100]); grid MINOR;
xlabel('Movement (s.d.)','FontSize',14)
ylabel('No of earthquakes','FontSize',14)
title('Theoretical Chi Distribution','FontSize',18); 
set(gca, 'fontsize', 14); 

% Printing figure to file
fig_filename = [simul_tag_col '_niter_' num2str(0) '.png'];
F    = getframe(FigH);
imwrite(F.cdata, fig_filename, 'png')
savefig(FigH,[fig_filename(1:end-4) '.fig'])
    

%% Other figures
for k = 1:n_iter
   
    FigH = figure('Position', get(0, 'Screensize'));
     
    subplot(1,2,1);
    scatter3(x_colap(k,:), y_colap(k,:), z_colap(k,:),'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
    %xlim([-4 4]); ylim([-4 4]); zlim([-4 4]);
    xlabel('Longitude (km)','FontSize',14)
    ylabel('Latitude (km)','FontSize',14)
    zlabel('depth (km)','FontSize',14)
    title(['Collapsed Earthquake Distribution (n\_iter = ', num2str(k),')'],'FontSize',18); 
    set(gca, 'fontsize', 14); view([-75.9 68.791015625]);
    
    grid MINOR;
    
    subplot(1,2,2)
    plot(x,chiscaled,'LineWidth',1.5);
    ylim([0 100]); 
    hold on;
    
    vv=nSD_moved(k,:); 
    vvv=vv(vv>0);
    [a,b] = hist(vvv,50); 
    hold on 
    stairs(b,a,'r','LineWidth',1.5)
    xlabel('Movement (s.d.)','FontSize',14)
    ylabel('No of earthquakes','FontSize',14)
    title(['Distribution of Hypocenter Movement (n\_iter = ', num2str(k),')'],'FontSize',18); 
    set(gca, 'fontsize', 14); grid MINOR;

%     hold on;
%     rr = sqrt(xs.^2 + ys.^2 + zs.^2);
%     [a,b] = hist(rr,50); 
%     hold on 
%     stairs(b,a,'r','LineWidth',1.5)    
   
    
%     pd = makedist('Chisquare');
%     pdf_normal = pdf(pd,vv);
%     plot(vv,pdf_normal*100,'go','LineWidth',2)
    
    % Printing figure to file
    fig_filename = [simul_tag_col '_niter_' num2str(k) '.png'];
    F    = getframe(FigH);
    imwrite(F.cdata, fig_filename, 'png')
    savefig(FigH,[fig_filename(1:end-4) '.fig'])
end



















% figure; subplot(1,2,1); plot3(x_colap, y_colap, z_colap,'bo'); view([0 90]); 
% xlim([-4 4]);
% ylim([-3 3]);
% grid MINOR;


% %% Plot original and collapsed hypocenters with the probability functions
% for k = 1:n_iter
%     
%     %if k==n_iter
%         % Wrte collapsed hypocenters to file
%         outfile = [simul_tag '_niter_' num2str(k) '.txt'];%'_uncert_' num2str(uncert)
%         fid=fopen(outfile,'w');
%         for kk=1:length(x_colap(k,:))
% 
%             fprintf(fid,'%12.5f %12.5f %12.5f\n',[x_colap(k,kk) y_colap(k,kk) z_colap(k,kk)]);
% 
%         end
%     %end
%     
%     % Extract distances moved by each hypocenter at the current iteration
%     % relative to its original location. Remove those that have not moved.
%     a = dist_to_ref_loc(k,:);
%     a(a<0.0001) = [];
%     
%     
%     % Figure 1
%     FigH = figure('Position', get(0, 'Screensize'));
% 
%     pos1 = [0.05 0.1 0.4 0.8];
%     ax1 = subplot('Position',pos1);
% 
%     scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
%     xlabel('Longitude (km)','FontSize',14)
%     ylabel('Latitude (km)','FontSize',14)
%     zlabel('depth (km)','FontSize',14)
%     title(['Earthquake Distribution of CSZ (' simul_tag,')'],'FontSize',18); 
%     set(gca, 'fontsize', 14);
%     view([-75.9 68.791015625]);
% 
% 
%     pos1 = [0.53 0.1 0.4 0.8];
%     ax2 = subplot('Position',pos1);
% 
%     scatter3(x_colap(k,:),y_colap(k,:),z_colap(k,:),'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 .75 .75])
%     xlabel('Longitude (km)','FontSize',14)
%     ylabel('Latitude (km)','FontSize',14)
%     zlabel('depth (km)','FontSize',14)
%     %title(['Collapsed Earthquake Distribution ( no of iteration = ', num2str(k),' )'],'FontSize',18); 
%     title(['Collapsed Earthquake Distribution (n\_iter = ', num2str(k),')'],'FontSize',18); 
%     set(gca, 'fontsize', 14);
%     view([-75.9 68.791015625]);
%     
%     Link = linkprop([ax1, ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
%     setappdata(gcf, 'StoreTheLink', Link);
% 
%     % Printing figure to file
%     fig_filename = [simul_tag '_e' num2str(uncert) 'km_niter_' num2str(k) '.png'];
%     F    = getframe(FigH);
%     imwrite(F.cdata, fig_filename, 'png')
%     savefig(FigH,[fig_filename(1:end-4) '.fig'])
%     
% 
%     % Figure 2    
%     Fig2 = figure('Position', get(0, 'Screensize'));
%     
%     %subplot(1,3,1)
%     len = 0.27; y_start = 0.25;
%     
%     pos1 = [0.05 y_start len 2*len];
%     subplot('Position',pos1)
%     pd1 = createFit_pdf(a);
%     title('Probability Density','FontSize',18); 
%     set(gca, 'fontsize', 14);
%       
%     %subplot(1,3,2)
%     pos2 = [0.365 y_start len 2*len];
%     subplot('Position',pos2)
%     pd2 = createFit_cdf(a);
%     title('Cummulative Probability Density','FontSize',18); 
%     set(gca, 'fontsize', 14);
%     
%     
%     %subplot(1,3,3)
%     pos2 = [0.7 y_start len 2*len];
%     subplot('Position',pos2)
%     pd3 = createFit_pf(a);
%     title(['Probability Plot (n\_iter = ', num2str(k),')'],'FontSize',18); 
%     set(gca, 'fontsize', 14);
%     
%     % Printing figure to file
%     fig_filename = [simul_tag '2_e' num2str(uncert) 'km_niter_' num2str(k) '.png'];
%     F2    = getframe(Fig2);
%     imwrite(F2.cdata, fig_filename, 'png')
%     savefig(Fig2,[fig_filename(1:end-4) '.fig'])
%     
%     a=[];
%          
% end
% 
