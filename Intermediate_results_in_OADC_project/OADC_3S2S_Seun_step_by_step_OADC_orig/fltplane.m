function [L,W,Strike,Dip,xv,yv,zv] = fltplane(V,D,xb,yb,zb)

% Calculate fault parameters from eigenvalues and eigenvectors of the
% Covariance matrix for this cluster

% eigenvalues and eigenvectors are sorted in ascending order

% L = fault length in the direction of max(D)
% W = fault width in the direction of intermediate(D)
% Strike, Dip, Rake of the plane
% xv(1:4), yv(1:4), zv(1:4) = vertices of the rectangular plane

% compute dip and strike
con=180.0./pi;
Dip=con.*acos(V(3,1));
Strike=con.*atan2(V(1,1),V(2,1));

% compute length and width
L=sqrt(12.*D(3,3));
W=sqrt(12.*D(2,2));

L2=L./2.0;
W2=W./2.0;

%L2=L;
%W2=W;

% compute vertice locations
xv(1)=W2.*V(1,2)+L2.*V(1,3) + xb;
yv(1)=W2.*V(2,2)+L2.*V(2,3) + yb;
zv(1)=W2.*V(3,2)+L2.*V(3,3) + zb;

xv(2)=W2.*V(1,2)-L2.*V(1,3) + xb;
yv(2)=W2.*V(2,2)-L2.*V(2,3) + yb;
zv(2)=W2.*V(3,2)-L2.*V(3,3) + zb;

xv(3)=-W2.*V(1,2)-L2.*V(1,3) + xb;
yv(3)=-W2.*V(2,2)-L2.*V(2,3) + yb;
zv(3)=-W2.*V(3,2)-L2.*V(3,3) + zb;

xv(4)=-W2.*V(1,2)+L2.*V(1,3) + xb;
yv(4)=-W2.*V(2,2)+L2.*V(2,3) + yb;
zv(4)=-W2.*V(3,2)+L2.*V(3,3) + zb;

return;
