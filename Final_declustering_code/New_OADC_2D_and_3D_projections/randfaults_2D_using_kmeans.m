function randfaults_2D_using_kmeans(n0,kthick)

%  Construct n0 randomn faults
%  just need the vertice locations

%  FAULT_FLAG = 0,  Use all hypocenters.  Initialization of the Process
%  FAULT_FLAG = kthick (from splitfault.m),  Split the thickest fault into
%               two random planes of length L/2

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

global xc yc vec_plane xb_old yb_old xs ys N Nc
global L xv yv L_old xv_old yv_old fscale
global xt yt Nt xb yb lambda3 Strike kmin

%fprintf('From randfaults\n');

xst=xt(kthick,1:Nt(kthick));
yst=yt(kthick,1:Nt(kthick));

X = [xst' yst'];

% 
% n0 = 2;
% X = [xs' ys'];

opts = statset('Display','final');
rng('default');
[idx,~] = kmeans(X,n0,'Distance','cityblock',... %'correlation' 'cityblock' 'sqeuclidean' 'cosine'
    'Replicates',5,'Options',opts);

%figure;
%colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'b.', 'k.', 'g.', 'c.', 'm.'};


    % Length, vec_plane, xv, yv
    for i=1:n0
        % compute the covariance matrix for this cluster
        xx = X(idx==i,1);
        yy = X(idx==i,2);
        
        xb(i) = mean(xx);
        yb(i) = mean(yy);
        
        Y = [xx yy];
        
        Cxy=cov(Y,0);

        % compute the eigenvalues and eigenvectors for this cluster
        [V,D]=eig(Cxy);

        % save the plane unit normal vector and eigenvalue
        vec_plane(i,1:2)=V(1:2,1);

        % calculate fault plane parameters from the eigen results
        % and calculate the vertices of the fault plane    
        [L(i),Strikek(i),xv(i,:),yv(i,:)] = fltplane_2D(Y,V,D,xb(i),yb(i));      

       % plot(xx,yy,colors{i}); hold on;
    end
    
    
    
%    
%     
%     %  plot final planes
% picname='Final Model';
% %datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
% datplot_2D(xs,ys,2,xv,yv,picname,'rr'); 

    
    

end




    
    
