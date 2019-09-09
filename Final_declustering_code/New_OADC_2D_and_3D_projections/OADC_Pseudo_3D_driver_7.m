%function OADC_Pseudo_3D_driver
%(simulation_tag)
%close all; clc; clear all;
tic

global orig_xs orig_ys orig_zs xs ys
global fscale N_thresh
global err_av kmin kmax N_loop simul_tag infile
global use_glo_var con_tol database database_lambda_only
global L_l  W_l  Strike_l  Dip_l  xv_l  yv_l  zv_l  vec_plane_l  lambda3_l
global L_ln W_ln Strike_ln Dip_ln xv_ln yv_ln zv_ln vec_plane_ln lambda3_ln
global lambda3 line_dens_incr

% ********************** Set Parameters ************************************
kmin = 1; kmax=20; err_av=1.5; %0.2 for synth
N_loop = 1; %simul_tag = char(simulation_tag); %
simul_tag = 'Simul.real.err1_5.incr20.no7'; 
use_glo_var = 1; N_thresh = 4;
infile = 'Simul.1_hypos.txt'; line_dens_incr = 2; theta_incr = 20;
%infile = 'testdata.txt';
%infile = 'cluster3.txt';
%infile = 'Simul.now_ALL_hypos_hypos.txt';
%infile = 'Simul.2_ALL_hypos_hypos.txt';
%infile = 'CSZ_hypos.txt';
%infile = 'COLCUM.20F_hypos.txt';



az_array = 0:theta_incr:179; 
el_array = -90:theta_incr:90;

%   Fault length scale for random faults. Will be between 0 and fscale in km
fscale=50.0;
%   Convergence tolerance value for the clustering algorithm in
%   'faultcluster'.  Represents the smallest change in global variance with
%   hypocenter clustering iteration.  The clustering process will stop once
%   the change in global variance with iteration drops to this value or
%   smaller.
con_tol=0.01;  %  units usually in km

%***************** Initial Checks ****************************
rng('shuffle');

% remove previous simulations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

% Save diary
diaryname = [simul_tag '.myDiaryFile.txt'];
diary(diaryname) 

%***************** Read Catalog of Hypocenters ****************************
read_catalog_P3D(infile,simul_tag,1);
      
add_array = ones(length(orig_xs),1);
database = [orig_xs' orig_ys' orig_zs' 100*add_array add_array];
database_lambda_only = [orig_xs' orig_ys' orig_zs' 100*add_array  add_array];

ncount=0; total_count = length(az_array)*length(el_array);

%textprogressbar('Determining the best fault model: '); 

for az = az_array
    for el = el_array
        ncount = ncount+1;
        
        %***************** Read Catalog of Hypocenters ********************
        read_catalog_P3D(infile,simul_tag,0);

        % project the 3D hypos on 2D planes specified by the azimuth and
        % elevations of viewpoint
        [xs,ys] = proj_of_3D_to_2D_plane(orig_xs,orig_ys,orig_zs,az,el);    
        
        % Peform OADC_2D on the projected hypocenters
        OADC_2D_on_proj_hypos() 
               
        % Classifying the cluster, and assigning lambda2 to each hypocenter 
        classifying_clusters_from_OADC_2D()
        
        % Progress...
        perc = (ncount/total_count)*100;
        %textprogressbar(perc);        
    end
end

%textprogressbar('done');

% fit_planes_and_plot_clusters
fit_planes_and_plot_clusters_based_on_lambda2_only()
fit_planes_and_plot_clusters_based_on_lambda2_and_Neqs()

simul_time = toc;

% saving all variables to file
savevar_filename = [simul_tag '.saved_variables.mat'];
save(savevar_filename)

% Move all figures to a folder name with the simul_tag
eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',simul_tag, '*',simul_tag,'_results'))
