function initial_faults_using_kmeans_2D(n0)

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

X = [xs' ys'];

clustering_method = 'dbscan'; % 'dbscan' 'kmeans'

if  strcmp(clustering_method,'dbscan') 
    rng('default');
    idx = rev_dbscan(X,0.7,10); % The default distance metric is Euclidean distance
    %idx = DBSCAN(X,1,10); % The default distance metric is Euclidean distance

    %gscatter(X(:,1),X(:,2),idx); 
    %title('DBSCAN Using Euclidean Distance Metric')

    figure;
    colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'b.', 'k.', 'g.', 'c.', 'm.', 'r.', 'b.', 'g.', 'k.', 'b.', 'k.', 'g.', 'c.', 'm.', };

    n0 = max(idx);
    kmin = max(idx);
    
    for i=1:n0
        plot(X(idx==i,1),X(idx==i,2),colors{mod(i,11)+1},'MarkerSize',12); hold on

        % Calculate new barycenters for each cluster
        xb(i) = mean(X(idx==i,1));
        yb(i) = mean(X(idx==i,2));
        Nt(i)= length(X(idx==i,1));

    end
    title('DBSCAN Using Euclidean Distance Metric')
    
    % Length, vec_plane, xv, yv
    for i=1:n0
       
        % compute the covariance matrix for this cluster
        Y = [X(idx==i,1) X(idx==i,2)];
        Cxy=cov(Y,0);

        % compute the eigenvalues and eigenvectors for this cluster
        [V,D]=eig(Cxy);

        % save the plane unit normal vector and eigenvalue
        vec_plane(i,1:2)=V(1:2,1);

        % calculate fault plane parameters from the eigen results
        % and calculate the vertices of the fault plane    
        [L(i),Strike(i),xv(i,:),yv(i,:)] = fltplane_2D(Y,V,D,xb(i),yb(i));      

    end


elseif strcmp(clustering_method,'kmeans') 

    opts = statset('Display','final');
    rng('default');
    [idx,C] = kmeans(X,n0,'Distance','cityblock',... %'correlation' 'cityblock' 'sqeuclidean'
        'Replicates',5,'Options',opts);

    figure;
    colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'b.', 'k.', 'g.', 'c.', 'm.'};

    for i=1:n0
        plot(X(idx==i,1),X(idx==i,2),colors{i},'MarkerSize',12); hold on
    end
    plot(C(:,1),C(:,2),'kx', 'MarkerSize',15,'LineWidth',3)
    title('Cluster Assignments and Centroids'); hold off

    % Calculate new barycenters for each cluster
    xb(1:n0) = C(:,1)';
    yb(1:n0) = C(:,2)';

    for i=1:n0
        Nt(i)= length(X(idx==i,1));
    end

    % Length, vec_plane, xv, yv
    for i=1:n0
        % compute the covariance matrix for this cluster
        Y = [X(idx==i,1) X(idx==i,2)];
        Cxy=cov(Y,0);

        % compute the eigenvalues and eigenvectors for this cluster
        [V,D]=eig(Cxy);

        % save the plane unit normal vector and eigenvalue
        vec_plane(i,1:2)=V(1:2,1);

        % calculate fault plane parameters from the eigen results
        % and calculate the vertices of the fault plane    
        [L(i),Strike(i),xv(i,:),yv(i,:)] = fltplane_2D(Y,V,D,xb(i),yb(i));      

    end

end

return;
    
    
