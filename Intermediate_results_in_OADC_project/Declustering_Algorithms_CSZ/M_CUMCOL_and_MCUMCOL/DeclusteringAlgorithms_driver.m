% Declustering Algorithms driver
clear all; close all; clc;

%% Decluster using cumulative distribution of tetrahedra volume
infile = 'CSZ_hypos.txt';
simul_tag = 'SCUM';

PROB = 0.05;
print_clus_t0_file = 1;

% remove previous calculations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

% run decluster code
decluster_w_cum_vol(infile, PROB, simul_tag, print_clus_t0_file)

%%
% % Method 2
infile = [simul_tag '_hypos.txt'];%'CSZ_hypos.txt';
simul_tag_col = 'SCUMCOL';

n_iter = 4; uncert = 1.5; % parameters for collapsing algorithm
skim = 'Nichol_etal'; % Choose the skim to use.
             % 'JS' = Jone and Stewart (1997) (Centroid = mean)
             % 'Nichol_etal' = Nicholson et al.(2000) (Centroid = Weighted mean)
           
decluster_w_collapsing(infile, simul_tag_col, n_iter, uncert, skim)

%% Move all figures to a folder name with the simul_tag
eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',simul_tag, '*',simul_tag,'_results'))
