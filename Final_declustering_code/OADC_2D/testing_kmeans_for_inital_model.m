close all; clc; clear all;
tic

global xs ys fscale use_glo_var con_tol zs
global FM_file dist2FM_threshold dip_threshold

kmin = 1; kmax=10; err_av=0.5;
N_loop = 15; simul_tag = 'Simul.2.OADC2D'; use_glo_var = 2;
FM_file='FM_dataset.csv'; dist2FM_threshold = 1; dip_threshold = 10;
%infile = 'proj.40.25.txt';
%infile = 'Simul.1_hypos.txt';
infile = 'CSZ_hypos.txt';


rng('shuffle');

% remove previous simulations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

% Save diary
diaryname = [simul_tag '.myDiaryFile.txt'];
diary(diaryname) 

%% ********************** Set Parameters ************************************
%   Fault length scale for random faults. Will be between 0 and fscale in km
fscale=20.0;
%   Convergence tolerance value for the clustering algorithm in
%   'faultcluster'.  Represents the smallest change in global variance with
%   hypocenter clustering iteration.  The clustering process will stop once
%   the change in global variance with iteration drops to this value or
%   smaller.
con_tol=0.01;  %  units usually in km
PLOT_FLAG1=1; % =0, no intermediate loop plots of data and planes

%***************** Read Catalog of Hypocenters ****************************
%read_catalog_2D(infile,simul_tag);

read_catalog(infile,simul_tag);





%X = [xs' ys'];
X = [xs' ys' zs'];


clustering_method = 'dbscan'; % 'dbscan' 'kmeans'




if  strcmp(clustering_method,'dbscan') 
    
    idx = rev_dbscan(X,1.5,10); % The default distance metric is Euclidean distance
    
    figure;
    colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'b.', 'k.', 'g.', 'c.', 'm.', 'r.', 'b.', 'g.', 'k.', 'b.', 'k.', 'g.', 'c.', 'm.', };

    n0 = max(idx);
    
    for i=1:n0
        plot3(X(idx==i,1),X(idx==i,2),X(idx==i,3),colors{i},'MarkerSize',12); hold on

    end  
        
        
        
    gscatter(X(:,1),X(:,2),X(:,3),idx); 
    title('DBSCAN Using Euclidean Distance Metric')
    shg
    
    
    
    
elseif strcmp(clustering_method,'kmeans') 
    
    
    
    
    
opts = statset('Display','final');
rng('default');
[idx,C] = kmeans(X,10,'Distance','cityblock',... %'correlation' 'cityblock' 'sqeuclidean'
    'Replicates',5,'Options',opts);

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12); hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12); hold on
plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',12); hold on
plot(X(idx==4,1),X(idx==4,2),'k.','MarkerSize',12); hold on
plot(X(idx==5,1),X(idx==5,2),'y.','MarkerSize',12); hold on
plot(X(idx==6,1),X(idx==6,2),'b.','MarkerSize',12); hold on
plot(X(idx==7,1),X(idx==7,2),'k.','MarkerSize',12); hold on
plot(X(idx==8,1),X(idx==8,2),'g.','MarkerSize',12); hold on
plot(X(idx==9,1),X(idx==9,2),'c.','MarkerSize',12); hold on
plot(X(idx==10,1),X(idx==10,2),'m.','MarkerSize',12); hold on;
plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',3)
title('Cluster Assignments and Centroids'); hold off

end


