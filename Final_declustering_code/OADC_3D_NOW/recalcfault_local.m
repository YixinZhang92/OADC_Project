function [Ltn,Wtn,Striket,Dipt,xvt,yvt,zvt,vec_plane_t] = recalcfault_local(xt,yt,zt,Nt,xbt,ybt,zbt,n0)
%  recalcfault - Compute the covariance matrix Cxy of each cluster, perform
%  principal components analysis, create best fault planes around each
%  cluster barycenter.

% n0 = number of fault planes
% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

% Seun changed the lambda3(k)=sqrt(12.*D(1,1)) to lambda3(k)=sqrt(D(1,1));
% Seun checks if Cxy contains NaN. This happens if no hypocenter is close to
% one of the splitted faults.
% Seun added hypos in cluster to the input parameters to fltplane.m

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global Strike Dip
%  Analyze each cluster
for k=1:n0
    %k
    % compute the covariance matrix for this cluster
    Cxy=cov( [xt(k,1:Nt(k))' yt(k,1:Nt(k))' zt(k,1:Nt(k))'],0);
    
    % Seun checks if Cxy contains NaN. This happens if no hypocenter is close to
    % one of the splitted faults.
    NrNaN = sum(isnan(Cxy(:)));
    if NrNaN > 0
        
        NrNaN;
        xt(k,1:Nt(k))';
        continue
    end
    
    % compute the eigenvalues and eigenvectors for this cluster
    [V,D]=eig(Cxy);
            
    % calculate fault plane parameters from the eigen results
    % and calculate the vertices of the fault plane
    X = [xt(k,1:Nt(k))' yt(k,1:Nt(k))' zt(k,1:Nt(k))'];
    
    [Ltn(k),Wtn(k),Striket(k),Dipt(k),xvt(k,:),yvt(k,:),zvt(k,:)] = fltplane(X,V,D,xbt(k),ybt(k),zbt(k));
            
    % save the plane unit normal vector and eigenvalue
    vec_plane_t(k,1:3)=V(1:3,1);
    %lambda3(k)=sqrt(12.*D(1,1));
    lambda3_t(k)=sqrt(D(1,1));
    
            
end
Strike;

end

