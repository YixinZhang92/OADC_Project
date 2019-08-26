function test_rectdist

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb
global L W xv yv zv Lold Wold xv_old yv_old zv_old

%***********  Initialize various matrices and parameters ******************
nmax=1;
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

%**************************************************************************

N=2;
xs=[2 3];
ys=[2 3];
zs=[2 3];

randfaults(1);

m=1;
for k=1:2;
dmin=rectdist(k,m);
end;

datplot(xs,ys,zs,1,xv,yv,zv);



%**************************************************************************