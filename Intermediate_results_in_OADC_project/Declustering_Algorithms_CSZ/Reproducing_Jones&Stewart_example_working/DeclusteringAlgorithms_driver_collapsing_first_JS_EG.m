% Declustering Algorithms driver
clear all; close all; clc;
% (Jones and Stewart Example)

%% Import hypocenters
% synthetic hypocenters
mu= 0; sigma=1; nhypos = 1000;
outfile='rnd_hypos_3D.txt';
[data,R] = create_syn_rnd_hypo_3D(nhypos,mu,sigma,outfile);







%% Create theoretical chi distribution
[x,chiscaled] = create_theor_chi_dist(nhypos,3);

% figures
figure;
subplot(1,2,1)
plot3(data(:,1),data(:,2),data(:,3),'ro'); grid MINOR;
view([0 90]); xlim([-3 3]); ylim([-3 3]);

subplot(1,2,2)
plot(x,chiscaled); shg
ylim([0 100]); grid MINOR
shg


%% Decluster using collapsing method first

n_iter = 5; nSD = 4; sigma_x = 1; sigma_y = 1; sigma_z = 1;% parameters for collapsing algorithm
skim = 'JS'; % Choose the skim to use.
             % 'JS' = Jone and Stewart (1997) (Centroid = mean)
             % 'Nichol_etal' = Nicholson et al.(2000) (Centroid = Weighted mean)
infile = outfile;
simul_tag_col = 'COL';

[x_colap, y_colap, z_colap, nSD_moved] = decluster_w_collapsing_JS_EG_3D...
(infile, simul_tag_col, n_iter, nSD, sigma_x, sigma_y, sigma_z, skim);

for k = [1 3 5] %1:n_iter
    
    figure; 
    subplot(1,2,1); plot3(x_colap(k,:), y_colap(k,:), z_colap(k,:),'bo'); view([0 90]); 
    xlim([-3 3]);
    ylim([-3 3]);
    grid MINOR;
    
    subplot(1,2,2)
    plot(x,chiscaled);
    ylim([0 100]); grid MINOR;
    hold on;
    
    xx = sort(nSD_moved(k,:));
    Y = chi2cdf(xx,3);
    nor = 0.9*sqrt(nhypos)/((max(Y))^2)*ones(1,length(Y));
    Yscaled = nor.*Y;
 
    
    plot(xx,Yscaled,'ro');
%     pd = makedist('Normal');
%     xx = sort(nSD_moved(k,:));
%     pdf_normal = pdf(pd,xx);
%     plot(xx,pdf_normal,'LineWidth',2);

end














% %% Decluster using collapsing method first
% infile = 'CSZ_hypos.txt';
% simul_tag_col = 'COL';
% 
% 
% 
% % remove previous calculations with the same simul_tag
% %eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag_col, '*', '*~'))
% 
% n_iter = 4; uncert = 1.5; % parameters for collapsing algorithm
% skim = 'Nichol_etal'; % Choose the skim to use.
%              % 'JS' = Jone and Stewart (1997) (Centroid = mean)
%              % 'Nichol_etal' = Nicholson et al.(2000) (Centroid = Weighted mean)
%            
% decluster_w_collapsing(infile, simul_tag_col, n_iter, uncert, skim)
% 
% %% Decluster using cumulative distribution of tetrahedra volume
% niter_to_decluser = 3;
% 
% infile = [simul_tag_col '_niter_' num2str(niter_to_decluser) '.txt'];%'CSZ_hypos.txt';
% simul_tag = 'COLCUM.2S';
% 
% PROB = 0.05;
% print_clus_t0_file = 1;
% 
% % run decluster code
% decluster_w_cum_vol(infile, PROB, simul_tag, print_clus_t0_file)
% 
% 
% %% Move all figures to a folder name with the simul_tag
% eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
% eval(sprintf('%s%s%s %s%s','! mv ',simul_tag_col, '*',simul_tag,'_results'))
