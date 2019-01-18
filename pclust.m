function pclust(kn,n0,nnext)

% Find n0 clusters in the data using k-means technique modified for
% distance measured to closest plane

% n0 = number of fault clusters to determine
% nnext = cluster choice to split into new barycenter and old plane
% kn = = iteration number

% Nc = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% vec_plane(1:Nc,1:3) = unit vector for the cluster plane
% xc,yc,zc = location of hypocenters in each cluster

% global variables definitions
global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale

if n0 == 1;     %  Initial call to pclust - form barycenter of entire data set
    
    Nt(1)=N;
    xt(1,1:N)=xs(1:N);
    yt(1,1:N)=ys(1:N);
    zt(1,1:N)=zs(1:N);
    
    xb(1)=mean(xt(1,1:N));
    yb(1)=mean(yt(1,1:N));
    zb(1)=mean(zt(1,1:N));
    
else;   %n0 > 1
    
    vec=vec_plane;  % initialize vec
    
    % load the new barycenters with most of the old before expanding to a
    % new barycenter
    
    %kk=0;
    %if kn == 1; kf=n0-1; else; kf=n0; end;
    %for k=1:kf;
    %    if k ~= nnext;
    %        kk=kk+1;
    %        xb(kk)=xb_old(k);
    %        yb(kk)=yb_old(k);
    %        zb(kk)=zb_old(k);
    %        vec(kk,1:3)=vec_plane(k,1:3);
    %        L(kk)=Lold(k);
    %        W(kk)=Wold(k);
    %        xv(kk,1:4)=xv_old(k,1:4);
    %        yv(kk,1:4)=yv_old(k,1:4);
    %        zv(kk,1:4)=zv_old(k,1:4));
            
    %    end;
    %end;
    
    %  Choose a random hypocenter from the "nnext" cluster to start a new
    %  point barycenter
    itest=randperm(Nc(nnext));
    ir=itest(1);
    
    %  Put this as the new barycenter
    %xb(n0)=xc(nnext,ir);
    %yb(n0)=yc(nnext,ir);
    %zb(n0)=zc(nnext,ir);
    
    %xb;
    %yb;
    %zb;
  
    
    % construct the pseudo-plane vector (unit vector) to give the correct
    % distance to the point barycenter (1,1,1).  The modulus of the dot
    % product of this vector with the (hypocenter - barycenter) vector will give
    % the distance to the barycenter.
    %vec(n0,1:3)=1.0;
    
    
    %vec;
    % Now go through the entire catalog of hypocenters to determine new
    % clusters.
    
    Nt(1:n0)=0;
    %N;
    for k=1:N;
        
        for m=1:n0;
            
            dst(m)=rectdist(k,m);
            
        end;
        
        dst;
        [val,index]=min(dst);
        
        Nt(index)=Nt(index) + 1;
        xt(index,Nt(index))=xs(k);
        yt(index,Nt(index))=ys(k);
        zt(index,Nt(index))=zs(k);
        
    end;
    
    % Calculate new barycenters for each cluster
    
    for kk=1:n0;
        nclus=Nt(kk);
 
        xb(kk)=mean(xt(kk,1:nclus));
        yb(kk)=mean(yt(kk,1:nclus));
        zb(kk)=mean(zt(kk,1:nclus));
    end;
    
    
    %  clusters have been produced, return for tests
    return;
    
    
end;
    
return;
