function [J]=Copy_of_pcluster_using_lamda3(n0)
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

% global variables definitions
global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global Strike Dip

% Now go through the entire catalog of hypocenters to determine new
% clusters.
    
J=0.0;
Nt(1:n0)=0;

 % Seun initialized these variables because it later becomes 
 % inconsistent with the total number of earthquakes in the catalog
xt = zeros(n0,N);
yt = zeros(n0,N);
zt = zeros(n0,N);


for k=1:N   %per hypocenter
        
    for m=1:n0  %per fault plane

        dst(m)=rectdist(k,m);
            
    end
        
    dst;
    
    %  find the closest fault plane
    [val,index]=min(dst);
        
    Nt(index)=Nt(index) + 1;
    xt(index,Nt(index))=xs(k);
    yt(index,Nt(index))=ys(k);
    zt(index,Nt(index))=zs(k);
    
    %accumulate the global variance
    J=J+val.*val;
        
end



%% Instead of calculating global variance, Seun determined the cluster lamda3
% lamda3 instead of global variance
J=J./round(N);




%......  Copy_of_recalcfault(n0) ............
% the result of interest is lamda3
Copy_of_recalcfault(n0)

% 
% Dip
% Strike

%% Calculate new barycenters for each cluster
    
for kk=1:n0
    nclus=Nt(kk);
 
    xb(kk)=mean(xt(kk,1:nclus));
    yb(kk)=mean(yt(kk,1:nclus));
    zb(kk)=mean(zt(kk,1:nclus));
end
    
    

% I am avoiding horizontal faults. Horizontal faults will e rotated by 90
% deg about the barycenter
index_of_horizontal_faults = find(Dip<20 & Dip~=0);
LHF = length(index_of_horizontal_faults);

if LHF > 0
    LHF;
   
    
    
    
    
   % [J]=Copy_of_pcluster_using_lamda3(n0)
    
end

%  clusters have been produced, return for tests

return;
 