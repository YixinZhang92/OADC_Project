% Clustering analysis driver
clear all; close all; clc;

global xs ys zs xv yv zv 

%% Generate synthetic datasets with hypocenters around faults
%  construct set of random hypocenters on the x-y plane with an additional
%  random error in z.  The x-y plane is specified by -L/2 < x < L/2 and
%  -W/2 < y W/2.  nhypos is the number of hypocenters desired.
%  zerr_av is the standard deviation of z error
%  strike, dip, and rake are in the geographical coordinate system
%  (degrees)
%  rt is a constant translation vector of the origin.
%  outfile is the name of the output file of hypocenter locations
%
% Syntax: rand_hypos(outfile,L,W,zerr_av,nhypos,strike,dip,rake,rt)

% Fault 1
rand_hypos('syn_hypo1.txt',10,5,0.1,100,45,70,0,[0 4 -5]);

% Fault 2
%rand_hypos('syn_hypo2.txt',6,4,0.1,100,135,45,0,[0 0 -5]);
rand_hypos('syn_hypo2.txt',6,4,0.1,100,45,70,0,[0 -2 -5]);

% Fault 3
rand_hypos('syn_hypo3.txt',10,5,0.1,100,45,70,0,[0 -4 -5]);

% Combining text files
system('rm syn_hypo_all.txt')
system('touch syn_hypo_all.txt')
system('cat syn_hypo1.txt >> syn_hypo_all.txt')
system('cat syn_hypo2.txt >> syn_hypo_all.txt')
system('cat syn_hypo3.txt >> syn_hypo_all.txt')
system('rm syn_hypo1.txt syn_hypo2.txt syn_hypo3.txt')

close all; 

%***************** Read Catalog of Hypocenters ****************************
%read_catalog('syn_hypo_all.txt');
read_catalog('testdata.txt');

%% Run the clustering analysis code and plot decustered datsets
strike_incr = 5; dip_incr= 10; mult_incr = 1; 
width= 1; % width of the depth slider when determining the number of EQs in each block.
min_eqs_for_a_cluster = 30;

[kmin,analy, value_counts] = clustering_analysis...
    (xs, ys, zs, strike_incr, dip_incr, width, mult_incr, min_eqs_for_a_cluster);

%% Extract fault geometries from the results of clustering analysis code
faults_from_clustering_analy(kmin,analy, value_counts);

%% Plot the datasets and faults
picname='Initial Model';
datplot(xs,ys,zs,kmin,xv,yv,zv,picname);

% END