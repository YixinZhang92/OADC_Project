function dmin=rectdist_2D(k,m)

global xc yc vec_plane xb_old yb_old xs ys N Nc
global xt yt Nt xb yb lambda3
global L xv yv Lold xv_old yv_old fscale

%fprintf('From rectdist %g %g\n',[k m]);

% k = index of hypocenter wanted
% m = index of fault plane wanted

% Choose hypocenter vector
u=[xs(k) ys(k);xs(k) ys(k)];
u1=[xs(k) ys(k)];

% Construct vertices vector
v=[xv(m,1:2)' yv(m,1:2)'];

% Construct distance vector from each vertex
d=u-v;
d1=[d(1,1) d(1,2)];
d2=[d(2,1) d(2,2)];

d1n=norm(d1);
d2n=norm(d2);
  
L1=L(m);

% Find all angles and perpendicular distances to edge
% 12 Face
%gam12=acos(d12./(d1n.*d2n))
if abs(d1n + d2n)<=L1 % hypocenter exactly at the center of the line.
    alph12=0;
    beta12=0;
    gam12=pi;
else
    [alph12,beta12,gam12]=trisol(d2n,d1n,L1,'r');
    a12=d1n.*d2n.*sin(gam12)./L1;
end

% Find the perpendicular distance to the infinite plane
vec(1:2)=vec_plane(m,1:2);
u2=u1-[xb(m) yb(m)];
dperp=abs(dot(u2,vec));

% Find the position field for the hypocenter relative to the plane and find
% the minimum distance.
p2=pi/2.0;
dmin=dperp;
if (alph12 > p2); dmin=d1n; end
if (beta12 > p2); dmin=d2n; end
    
%dmin;
return;
