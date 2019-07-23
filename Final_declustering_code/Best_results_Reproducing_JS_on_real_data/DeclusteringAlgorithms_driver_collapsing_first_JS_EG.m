% Declustering Algorithms driver
clear all; close all; clc;

% parameters for collapsing algorithm
%n_iter = 5; nSD = 4; sigma_x = 1; sigma_y = 1; sigma_z = 1;
n_iter = 20; nSD = 4; sigma_x = 0.15; sigma_y = 0.15; sigma_z = 0.35;
%n_iter = 20; nSD = 4; sigma_x = 1; sigma_y = 1; sigma_z = 1;

skim = 'JS'; % Choose the skim to use.
             % 'JS' = Jone and Stewart (1997) (Centroid = mean)
             % 'Nichol_etal' = Nicholson et al.(2000) (Centroid = Weighted mean)

simul_tag_col = 'COLs1_5';

% remove previous calculations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag_col, '*', '*~'))

%% Import hypocenters
syn = 0; % syn = 0 - Real hypocenters
         %     = 1 - Synthetic hypocenters
if syn == 1
    % synthetic hypocenters
    mu= 0; sigma=1; nhypos = 1000;
    outfile='rnd_hypos_3D.txt';
    infile = outfile;
    [data,R] = create_syn_rnd_hypo_3D(nhypos,mu,sigma,outfile);

else
    % Real data
    infile = 'CSZ_hypos.txt';
   
    fid=fopen(infile,'r');
    [data,~]=fscanf(fid,'%g %g %g',[3,inf]);
    fclose(fid);

    data=data';
    
%     catalognotinPL17y8899 = readtable('catalog_not_in_PL17_y88_99.txt');
%     xs_8899 = table2array(catalognotinPL17y8899(:,4)); xs_8899n = (xs_8899 + 70.6)*75.778;
%     ys_8899 = table2array(catalognotinPL17y8899(:,3)); ys_8899n = (ys_8899 - 47.2)*111.1743;
%     zs_8899 = table2array(catalognotinPL17y8899(:,5));
%     data_8899 = [xs_8899n ys_8899n -zs_8899];
%     
%     data = [data; data_8899];
%         
%     catalognotinPL17y0017 = readtable('catalog_not_in_PL17_y00_17.txt');
%     xs_0017 = table2array(catalognotinPL17y0017(:,6)); xs_0017n = (xs_0017 + 70.6)*75.778;
%     ys_0017 = table2array(catalognotinPL17y0017(:,5)); ys_0017n = (ys_0017 - 47.2)*111.1743;
%     zs_0017 = table2array(catalognotinPL17y0017(:,7));
%     data_0017 = [xs_0017n ys_0017n -zs_0017];
%     
%     data = [data; data_0017];
    
    
    
    
    [nhypos,~]=size(data);
    
    
    % write synthetic data to outfile
    data_file = 'all_hypos.txt';
    fid=fopen(data_file,'w');
    for kk=1:nhypos

        fprintf(fid,'%12.5f %12.5f %12.5f\n',[data(kk,1) data(kk,2) data(kk,3)]);

    end

    infile = data_file;
end

%% Decluster using collapsing method first

[x_colap, y_colap, z_colap, nSD_moved] = decluster_w_collapsing_JS_EG_3D...
(infile, simul_tag_col, n_iter, nSD, sigma_x, sigma_y, sigma_z, skim);


%% Decluster using cumulative distribution of tetrahedra volume
niter_to_decluser = 20;

infile = [simul_tag_col '_niter_' num2str(niter_to_decluser) '.txt'];%'CSZ_hypos.txt';
simul_tag = 'COLCUM.20F';

PROB = 0.05;
print_clus_t0_file = 1;

% run decluster code
[vol_sort, volrand_sort,V05, NV05_cat] = decluster_w_cum_vol(infile, PROB, simul_tag, print_clus_t0_file);


%% Move all figures to a folder name with the simul_tag
eval(sprintf('%s%s%s','! mkdir ',simul_tag_col,'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',simul_tag_col, '*',simul_tag_col,'_results'))





























%     xx = sort(nSD_moved(k,:));
%     Y = chi2cdf(xx,3);
%     nor = 0.9*sqrt(nhypos)/((max(Y))^2)*ones(1,length(Y));
%     Yscaled = nor.*Y;
%     
%     plot(xx,Yscaled,'ro');
    
    
%     pd = makedist('Normal');
%     xx = sort(nSD_moved(k,:));
%     pdf_normal = pdf(pd,xx);
%     plot(xx,pdf_normal,'LineWidth',2);


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
