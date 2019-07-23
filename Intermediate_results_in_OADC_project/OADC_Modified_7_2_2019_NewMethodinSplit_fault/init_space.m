function init_space(kmax)
%  init_space - initialize working arrays for fault clusters
% Choice of cluster to make a new barycenter

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt_old yt_old zt_old vec_plane_old lambda3_old
global Strike Dip Strike_old Dip_old Nt_old
% cluster location matrices
xc(1:kmax,1:N)=0.;
yc(1:kmax,1:N)=0.;
zc(1:kmax,1:N)=0.;

% number of events in each cluster
Nc(1:kmax)=0;

% trial cluster location matrices
xt(1:kmax,1:N)=0.;
yt(1:kmax,1:N)=0.;
zt(1:kmax,1:N)=0.;

xt_old(1:kmax,1:N)=0.;
yt_old(1:kmax,1:N)=0.;
zt_old(1:kmax,1:N)=0.;

% trial cluster barycenter location matrices
xb(1:kmax)=0.;
yb(1:kmax)=0.;
zb(1:kmax)=0.;

xb_old(1:kmax)=0.;
yb_old(1:kmax)=0.;
zb_old(1:kmax)=0.;

% number of events in each trial cluster
Nt(1:kmax)=0;

Nt_old(1:kmax)=0;

% cluster covariance matrix
Cxy(1:kmax,1:3,1:3)=0.;

% eigenvector that describes each plane
vec_plane(1:kmax,1:3)=0.;
vec_plane_old(1:kmax,1:3)=0.;

% fault plane parameters
L(1:kmax)=0.;
W(1:kmax)=0.;
Strike(1:kmax)=0.;
Dip(1:kmax)=0.;

L_old(1:kmax)=0.;
W_old(1:kmax)=0.;
Strike_old(1:kmax)=0.;
Dip_old(1:kmax)=0.;

% minimum eigenvalue
lambda3(1:kmax)=0.;

lambda3_old(1:kmax)=0.;

% fault plane vertices for plot purposes
xv(1:kmax,1:4)=0.;
yv(1:kmax,1:4)=0.;
zv(1:kmax,1:4)=0.;

xv_old(1:kmax,1:4)=0.;
yv_old(1:kmax,1:4)=0.;
zv_old(1:kmax,1:4)=0.;

end

