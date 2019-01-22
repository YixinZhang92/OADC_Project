function faults_from_clustering_analy(n0,analy, value_counts)

%  Construct n0 faults based on the clustering method
%  just need the vertice locations

%  FAULT_FLAG = 0,  Use all hypocenters.  Initialization of the Process
%  FAULT_FLAG = kthick (from splitfault.m),  Split the thickest fault into
%               two random planes of length L/2

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster



global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt yt zt Nt xb yb zb lambda3

%fprintf('From randfaults\n');

for k=1:n0
    % extracting locations of each eq in this cluster
    x_clus = analy(analy(:,8) == value_counts(k,1),1);
    y_clus = analy(analy(:,8) == value_counts(k,1),2);
    z_clus = analy(analy(:,8) == value_counts(k,1),3);
    
    % trial cluster barycenter location matrices
    xb(k)=mean(x_clus);
    yb(k)=mean(y_clus);
    zb(k)=mean(z_clus);

    % compute the covariance matrix for this cluster
    Cxy=cov( [x_clus y_clus z_clus],0);
        
    % compute the eigenvalues and eigenvectors for this cluster
    [V,D]=eig(Cxy);
            
    % calculate fault plane parameters from the eigen results
    % and calculate the vertices of the fault plane
    [L(k),W(k),Strike(k),Dip(k),xv(k,:),yv(k,:),zv(k,:)] = fltplane(V,D,xb(k),yb(k),zb(k));
            
    % save the plane unit normal vector and eigenvalue
    vec_plane(k,1:3)=V(1:3,1);
    lambda3(k)=sqrt(12.*D(1,1));
    
end

return;
    
    
