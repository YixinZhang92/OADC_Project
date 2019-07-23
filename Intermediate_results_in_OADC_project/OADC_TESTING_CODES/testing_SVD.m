clear all; close all; clc;

global xs ys zs xv yv zv 

system('rm syn_hypo_all.txt')
system('touch syn_hypo_all.txt')

% Fault 1
%rand_hypos('syn_hypo1.txt',10,5,0.1,100,45,70,0,[0 4 -5]);
rand_hypos('syn_hypo1.txt',10,5,0.0,100,0,90,0,[0 4 -5]);
system('cat syn_hypo1.txt >> syn_hypo_all.txt')

% % % Fault 2
% % %rand_hypos('syn_hypo2.txt',6,4,0.1,100,135,45,0,[0 0 -5]);
% rand_hypos('syn_hypo2.txt',6,4,0.1,100,0,90,0,[2 4 -5]);
% system('cat syn_hypo2.txt >> syn_hypo_all.txt')
% % 
% % % Fault 3
% rand_hypos('syn_hypo3.txt',10,5,0.1,100,45,70,0,[0 -4 -5]);
% system('cat syn_hypo3.txt >> syn_hypo_all.txt')

% Combining text files
system('rm syn_hypo1.txt syn_hypo2.txt syn_hypo3.txt')


close all; 
%***************** Read Catalog of Hypocenters ****************************
%read_catalog('syn_hypo_all.txt');
read_catalog('syn_hypo_all.txt');


data = [xs' ys' zs'];

% [U,S,V] = svd(data);
k=1;
% trial cluster barycenter location matrices
xb=mean(xs');
yb=mean(ys');
zb=mean(zs');

% compute the covariance matrix for this cluster
Cxy=cov( [xs' ys' zs'],0);

% compute the eigenvalues and eigenvectors for this cluster
[V,D]=eig(Cxy);

% calculate fault plane parameters from the eigen results
% and calculate the vertices of the fault plane
[L(k),W(k),Strike(k),Dip(k),xv(k,:),yv(k,:),zv(k,:)] = fltplane(V,D,xb(k),yb(k),zb(k));

% save the plane unit normal vector and eigenvalue
vec_plane(k,1:3)=V(1:3,1);
lambda3(k)=sqrt(12.*D(1,1));