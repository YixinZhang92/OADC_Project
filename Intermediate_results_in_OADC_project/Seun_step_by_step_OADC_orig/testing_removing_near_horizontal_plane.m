

close all; clear all; clc;

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global Strike Dip

close all;
%kmin = 1;kmax=3;err_av=0.5;infile='testdata.txt';

kmin=1;kmax=1;err_av=0.05;infile= 'result_declustered_collapsed_hypo_0.8.txt';

%********************** Set Parameters ************************************
%   Fault length scale for random faults.  Will be between 0 and fscale in
%   km
fscale=50.0;
%   Convergence tolerance value for the clustering algorithm in
%   'faultcluster'.  Represents the smallest change in global variance with
%   hypocenter clustering iteration.  The clustering process will stop once
%   the change in global variance with iteration drops to this value or
%   smaller.
con_tol=0.001;  %  units usually in km
PLOT_FLAG1=1;   % =0, no intermediate loop plots of data and planes

%***************** Read Catalog of Hypocenters ****************************
read_catalog(infile);

%********************** Initialize Space **********************************
init_space(kmax);

%******************* Initialize random faults *****************************
FAULT_FLAG=0;   % Initialization, use all hypocenters
randfaults(kmin,FAULT_FLAG);

%  plot initial planes
picname='Initial Model';
datplot(xs,ys,zs,kmin,xv,yv,zv,picname);

xb
yb
zb



%%

xvn = xv-xb;
yvn = yv-yb;
zvn = zv-zb;

strike = 0;
dip = 90;

%% concatenating x,y and z
    R=[xvn;  yvn; zvn;];
    
    Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];
    Rdip= Ddip*R;

    % add bac the barycenter
    Rdip(1,:) = Rdip(1,:) + xb;
    Rdip(2,:) = Rdip(2,:) + yb;
    Rdip(3,:) = Rdip(3,:) + zb;
    
    xvn = Rdip(1,:);
    yvn = Rdip(2,:);
    zvn = Rdip(3,:);
    
  
    datplot(xs,ys,zs,kmin,xvn,yvn,zvn,picname);