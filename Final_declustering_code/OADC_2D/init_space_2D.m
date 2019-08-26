function init_space_2D(kmax)
%  init_space - initialize working arrays for fault clusters
% Choice of cluster to make a new barycenter

global xc yc vec_plane xb_old yb_old xs ys N Nc
global xt yt Nt xb yb lambda3
global L xv yv L_old xv_old yv_old fscale
global xt_old yt_old vec_plane_old lambda3_old
global Strike Strike_old Nt_old
% cluster location matrices
xc(1:kmax,1:N)=0.;
yc(1:kmax,1:N)=0.;

% number of events in each cluster
Nc(1:kmax)=0;

% trial cluster location matrices
xt(1:kmax,1:N)=0.;
yt(1:kmax,1:N)=0.;

xt_old(1:kmax,1:N)=0.;
yt_old(1:kmax,1:N)=0.;

% trial cluster barycenter location matrices
xb(1:kmax)=0.;
yb(1:kmax)=0.;

xb_old(1:kmax)=0.;
yb_old(1:kmax)=0.;

% number of events in each trial cluster
Nt(1:kmax)=0;

Nt_old(1:kmax)=0;

% cluster covariance matrix
Cxy(1:kmax,1:2,1:2)=0.;

% eigenvector that describes each plane
vec_plane(1:kmax,1:2)=0.;
vec_plane_old(1:kmax,1:2)=0.;

% fault plane parameters
L(1:kmax)=0.;
Strike(1:kmax)=0.;

L_old(1:kmax)=0.;
Strike_old(1:kmax)=0.;

% minimum eigenvalue
lambda3(1:kmax)=0.;

lambda3_old(1:kmax)=0.;

% fault plane vertices for plot purposes
xv(1:kmax,1:2)=0.;
yv(1:kmax,1:2)=0.;

xv_old(1:kmax,1:2)=0.;
yv_old(1:kmax,1:2)=0.;

end

