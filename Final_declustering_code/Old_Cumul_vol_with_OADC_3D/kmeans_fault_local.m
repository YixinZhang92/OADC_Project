function [xbt,ybt,zbt,xvt,yvt,zvt,vec_plane_t,Ltn,Wtn] = kmeans_fault_local(xst,yst,zst,n0)

X = [xst' yst' zst'];

opts = statset('Display','final');
[idx,C] = kmeans(X,2,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

xbt = C(:,1);
ybt = C(:,2);
zbt = C(:,3);

for k=1:n0
     % compute the covariance matrix for this cluster
    Cxy=cov(X(idx==k,:),0);
    
    % Seun checks if Cxy contains NaN. This happens if no hypocenter is close to
    % one of the splitted faults.
    NrNaN = sum(isnan(Cxy(:)));
    if NrNaN > 0
        continue
    end
    
    % compute the eigenvalues and eigenvectors for this cluster
    [V,D]=eig(Cxy);
            
    % calculate fault plane parameters from the eigen results
    % and calculate the vertices of the fault plane
    XX = X(idx==k,:);
    
    [Ltn(k),Wtn(k),Striket(k),Dipt(k),xvt(k,:),yvt(k,:),zvt(k,:)] = fltplane(XX,V,D,xbt(k),ybt(k),zbt(k));
            
    % save the plane unit normal vector and eigenvalue
    vec_plane_t(k,1:3)=V(1:3,1);
    %lambda3(k)=sqrt(12.*D(1,1));
    lambda3_t(k)=sqrt(D(1,1));
         
end  
    
return;
    
    
% figure;
% plot3(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',12)
% hold on
% plot3(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',12)
% plot3(C(:,1),C(:,2),C(:,3),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% title 'Cluster Assignments and Centroids'
% hold off