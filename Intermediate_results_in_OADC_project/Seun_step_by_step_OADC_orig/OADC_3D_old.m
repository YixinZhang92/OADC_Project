function OADC_3D_old(nmax,err_av,infile)

%  Implementation of 3-D Optimal Anisotropic Dynamic Clustering from

%  Ouillon, Ducorbier, and Sornette (2008).  Automatic reconstruction of
%  fault networks from seismicity catalogs: Three-dimensional optimal
%  anisotropic dynamic clustering, JGR, 113, B01306,
%  doi:10.1029/2007JB00503.

%  specify:

%       nmax = maximum number of fault planes analyzed
%       err_av = average hypocentral error in km
%       infile = file containing (x,y,z) positions of hypocenters

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb
global L W xv yv zv Lold Wold xv_old yv_old zv_old fscale

close all;

%****************  Input hypocentral locations ****************************
fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

% plot the input data
figure;
plot3(xs,ys,zs,'o');
axis equal;
title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

%***********  Initialize various matrices and parameters ******************

% Choice of cluster to make a new barycenter
nnext=0;

% cluster location matrices
xc(1:nmax,1:N)=0.;
yc(1:nmax,1:N)=0.;
zc(1:nmax,1:N)=0.;

% number of events in each cluster
Nc(1:nmax)=0;

% trial cluster location matrices
xt(1:nmax,1:N)=0.;
yt(1:nmax,1:N)=0.;
zt(1:nmax,1:N)=0.;

% trial cluster barycenter location matrices
xb(1:nmax)=0.;
yb(1:nmax)=0.;
zb(1:nmax)=0.;

xb_old(1:nmax)=0.;
yb_old(1:nmax)=0.;
zb_old(1:nmax)=0.;

% number of events in each trial cluster
Nt(1:nmax)=0;

% cluster covariance matrix
Cxy(1:nmax,1:3,1:3)=0.;

% eigenvector that describes each plane
vec_plane(1:nmax,1:3)=0.

% fault plane parameters
L(1:nmax)=0.;
W(1:nmax)=0.;
Strike(1:nmax)=0.;
Dip(1:nmax)=0.;
Rake(1:nmax)=0.;

% minimum eigenvalue
lambda3(1:nmax)=0.;

% fault plane vertices for plot purposes
xv(1:nmax,1:4)=0.;
yv(1:nmax,1:4)=0.;
zv(1:nmax,1:4)=0.;

% scale of fault plane length used in assigning initial random faults to
% the catalog (km)
fscale=10.0;

%**************************************************************************

%  Now for the big loop over the number of possible fault planes considered
%  (read n-naught or n-zero, n0 is the number of planes considered in each 
%   iteration)

nfinish=0;

for n0=1:nmax;
    
    fprintf('Planes %d\n',n0);
    
    %  Assume niter # of iterations when n0 > 1
    if n0 == 1; niter=1; else; niter=40; end;
    
    for kk=1:niter;
        
        fprintf('Iteration %d\n',kk);
        
    if nfinish == 0
        
        %  construct n0 random fault planes for cluster seeds if kk=1 to 
        %  start off interations.  Will compute L, W, and vertices vectors
            if kk == 1 && n0 > 1;
                randfaults(n0);
            end;
            
            %  Find new clusters
            pclust(kk,n0,nnext);
            
    
        %  Analyze each cluster
        for k=1:n0;
            %k
            % compute the covariance matrix for this cluster
            Cxy=cov( [xt(k,1:Nt(k))' yt(k,1:Nt(k))' zt(k,1:Nt(k))'],0);
        
            % compute the eigenvalues and eigenvectors for this cluster
            [V,D]=eig(Cxy);
            
            % calculate fault plane parameters from the eigen results
            % and calculate the vertices of the fault plane
            [L(k),W(k),Strike(k),Dip(k),xv(k,:),yv(k,:),zv(k,:)] = fltplane(V,D,xb(k),yb(k),zb(k));
            
            % save the plane unit normal vector and eigenvalue
            vec_plane(k,1:3)=V(1:3,1);
            lambda3(k)=sqrt(12.*D(1,1));
            
        end;
    
        %  plot results for this iteration
        if nfinish == 1 || kk == niter;
            datplot(xs,ys,zs,n0,xv,yv,zv);
        end;
        
        %  replace old clusters with new clusters
        xc=xt;
        yc=yt;
        zc=zt;
        Nc=Nt;
        xb_old=xb;
        yb_old=yb;
        zb_old=zb;
        
        %  Test all planes for thickness, choose maximum thickness for next
        %  iteration
        [thick,nnext]=max(lambda3);
    
        %  print results to the screen
        fprintf('Thickness of each plane\n');
        for k=1:n0; fprintf('%g %g %g\n',lambda3(k));end;
    
        if thick <= err_av;
        
            nfinish=1;
            fprintf('Clustering Complete - Ending\n');
        end;
    
    else
    end;    %nfinish loop
    
    end;    %niter loop
        
end;    %n0 loop
