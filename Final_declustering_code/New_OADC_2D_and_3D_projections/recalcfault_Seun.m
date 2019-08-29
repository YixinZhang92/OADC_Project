function recalcfault_Seun(n0)
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
        
    % Seun changed the way we determine the fault vertex, length and width
    % of the fault.
    % calculate fault plane parameters from the eigen results
    % and calculate the vertices of the fault plane
    
    % compute vertice locations
    xv(k,1)=min(xt(k,1:Nt(k)));
    yv(k,1)=max(yt(k,1:Nt(k)));

    xv(k,2)=min(xt(k,1:Nt(k)));
    yv(k,2)=min(yt(k,1:Nt(k)));

    xv(k,3)=max(xt(k,1:Nt(k)));
    yv(k,3)=min(yt(k,1:Nt(k)));

    xv(k,4)=max(xt(k,1:Nt(k)));
    yv(k,4)=max(yt(k,1:Nt(k)));
    
    zv(k,:) = (1/V(3,1)) .* ((xb(k)*V(1,1) + yb(k)*V(2,1) + (zb(k)*V(3,1))) ...
        - (V(1,1)*[xv(k,1) xv(k,2) xv(k,3) xv(k,4)] + V(2,1)*[yv(k,1) yv(k,2) yv(k,3) yv(k,4)] ));

    L(k) = sqrt((xv(k,2)-xv(k,1))^2 + (yv(k,2)-yv(k,1))^2 + (zv(k,2)-zv(k,1))^2);
    W(k) = sqrt((xv(k,4)-xv(k,1))^2 + (yv(k,4)-yv(k,1))^2 + (zv(k,4)-zv(k,1))^2);
    
    % compute dip and strike
    con=180.0./pi;
    % Dip=con.*acos(V(3,1));
    % Strike=con.*atan2(V(1,1),V(2,1));

    Dip(k)=con.*acos(V(3,1));
    if (V(1,1) >= 0) && (V(2,1) >= 0)
        Strike(k)=con.*atan(abs(V(1,1))/abs(V(2,1)))-90;
    elseif (V(1,1) >= 0) && (V(2,1) < 0)
        Strike(k)=180-con.*atan(abs(V(1,1))/abs(V(2,1)))-90;
    elseif (V(1,1) < 0) && (V(2,1) < 0)
        Strike(k)=180+con.*atan(abs(V(1,1))/abs(V(2,1)))-90;    
    elseif (V(1,1) < 0) && (V(2,1) >= 0)
        Strike(k)=360-con.*atan(abs(V(1,1))/abs(V(2,1)))-90;    
    end    

    if Strike(k) < 0
        Strike(k) = Strike(k)+180;
    end
      
    % save the plane unit normal vector and eigenvalue
    vec_plane(k,1:3)=V(1:3,1);
    %lambda3(k)=sqrt(12.*D(1,1));
    lambda3(k)=sqrt(D(1,1));
              
end
Strike;

end

