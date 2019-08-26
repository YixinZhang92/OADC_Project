function recalcfault_2D(n0)
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

global xc yc vec_plane xb_old yb_old xs ys N Nc
global xt yt Nt xb yb lambda3
global L xv yv L_old xv_old yv_old fscale
global Strike
%  Analyze each cluster
for k=1:n0
    %k
    % compute the covariance matrix for this cluster
    Cxy=cov( [xt(k,1:Nt(k))' yt(k,1:Nt(k))'],0);
    
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
    X = [xt(k,1:Nt(k))' yt(k,1:Nt(k))'];
    
    [L(k),Strike(k),xv(k,:),yv(k,:)] = fltplane_2D(X,V,D,xb(k),yb(k));
            
    % save the plane unit normal vector and eigenvalue
    vec_plane(k,1:2)=V(1:2,1);
    %lambda3(k)=sqrt(12.*D(1,1));
    lambda3(k)=sqrt(D(1,1));
    
    
    % plot density 
    
    [nhypos,~] = size(X);
    projX=[V(:,2) V(:,1)]\X';
    pxs = projX(1,:);
    pys = projX(2,:);
    
    figure;
    hist(pxs, min(pxs):1:max(pxs));shg
    
    
%         for i = 1:nhypos
%     AP = [0 0] + X(i,:);
%     AB = [0 0] - [V(1,2) V(2,2)];
%     NN(i,:) = [0 0] + dot(AP,AB) / dot(AB,AB) * AB;
% 
%         end
%     
%         

  
end
Strike;

end

