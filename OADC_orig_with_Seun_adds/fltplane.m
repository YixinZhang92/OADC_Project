function [L,W,Strike,Dip,xv,yv,zv] = fltplane(V,D,xb,yb,zb)

% Calculate fault parameters from eigenvalues and eigenvectors of the
% Covariance matrix for this cluster

% eigenvalues and eigenvectors are sorted in ascending order

% L = fault length in the direction of max(D)
% W = fault width in the direction of intermediate(D)
% Strike, Dip, Rake of the plane
% xv(1:4), yv(1:4), zv(1:4) = vertices of the rectangular plane

% Seun changed the equations for determining strike and dip accordingly for
% strike to go from the north clockwisely.

% compute dip and strike
con=180.0./pi;
% Dip=con.*acos(V(3,1));
% Strike=con.*atan2(V(1,1),V(2,1));

Dip=con.*acos(V(3,1));
if (V(1,1) >= 0) && (V(2,1) >= 0)
    Strike=con.*atan(abs(V(1,1))/abs(V(2,1)))-90;
elseif (V(1,1) >= 0) && (V(2,1) < 0)
    Strike=180-con.*atan(abs(V(1,1))/abs(V(2,1)))-90;
elseif (V(1,1) < 0) && (V(2,1) < 0)
    Strike=180+con.*atan(abs(V(1,1))/abs(V(2,1)))-90;    
elseif (V(1,1) < 0) && (V(2,1) >= 0)
    Strike=360-con.*atan(abs(V(1,1))/abs(V(2,1)))-90;    
end    
    
if Strike < 0
    Strike = Strike+180;
end

% if V(2,1) >= 0
%     Strike=con.*atan2(V(1,1),V(2,1))+90;
% else
%     Strike=con.*atan2(V(1,1),V(2,1))+270;
% end

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
