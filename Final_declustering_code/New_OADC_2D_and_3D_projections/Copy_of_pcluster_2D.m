function [J]=pcluster(n0)
%  pcluster - form n0 clusters of seismicity around the input fault planes
%  using the modified k-means method.
%  Also compute the global variance of the fit.

% n0 = number of fault clusters to determine
% J  = global variance is output

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

% Seun initialized xt, yt and zt variables because it later becomes 
% inconsistent with the total number of earthquakes in the catalog
% Moved recalcfault.m from faultcluster.m to Copy_of_pcluster.m


% global variables definitions
global xc yc vec_plane xb_old yb_old xs ys N Nc
global xt yt Nt xb yb lambda3
global L xv yv L_old xv_old yv_old fscale

% Now go through the entire catalog of hypocenters to determine new
% clusters.
    
J=0.0;
Nt(1:n0)=0;
%N;
% Seun initialized these variables because it later becomes 
% inconsistent with the total number of earthquakes in the catalog
xt = zeros(n0,N);
yt = zeros(n0,N);

for k=1:N   %per hypocenter
    for m=1:n0  %per fault plane
        dst(m)=rectdist_2D(k,m);         
    end
        
    dst;
    
    %  find the closest fault plane
    [val,index]=min(dst);
        
    Nt(index)=Nt(index) + 1;
    xt(index,Nt(index))=xs(k);
    yt(index,Nt(index))=ys(k);
    
    %  accumulate the global variance
    J=J+val.*val;
        
end

%  global variance
J=J./round(N);

% Seun moved this here.
%  Compute Cxy matrix, perform principal components analysis, create
%  new fault planes for each cluster

if min(Nt(1:n0)) > 1
    
    recalcfault_2D(n0);
else
    
    J = 10000; % make it big so that the model could be rejected
end

% Calculate new barycenters for each cluster
for kk=1:n0
    nclus=Nt(kk);

    xb(kk)=mean(xt(kk,1:nclus));
    yb(kk)=mean(yt(kk,1:nclus));
   
end
       
%  clusters have been produced, return for tests
return;
 