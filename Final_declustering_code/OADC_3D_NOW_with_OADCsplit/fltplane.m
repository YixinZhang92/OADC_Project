function [L,W,Strike,Dip,xv,yv,zv] = fltplane(X,V,D,xb,yb,zb)

% Calculate fault parameters from eigenvalues and eigenvectors of the
% Covariance matrix for this cluster

% eigenvalues and eigenvectors are sorted in ascending order

% L = fault length in the direction of max(D)
% W = fault width in the direction of intermediate(D)
% Strike, Dip, Rake of the plane
% xv(1:4), yv(1:4), zv(1:4) = vertices of the rectangular plane

% Seun changed the equations for determining strike and dip accordingly for
% strike to go from the north clockwisely.
% Seun added hypos in cluster to the input parameters to fltplane.m
% Seun changed the way to determine L and W of fault plane

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
use_L_Ouillon = 1;

if use_L_Ouillon == 1
    
    L=sqrt(12.*D(3,3));
    W=sqrt(12.*D(2,2));

else

    [nhypos,~] = size(X);

    for i = 1:nhypos
    AP = [0 0 0] + X(i,:);
    AB = [0 0 0] - [V(1,3) V(2,3) V(3,3)];
    NN(i,:) = [0 0 0] + dot(AP,AB) / dot(AB,AB) * AB;

    end
    L = sqrt((max(NN(:,1))-min(NN(:,1)))^2 + (max(NN(:,2))-min(NN(:,2)))^2 +(max(NN(:,3))-min(NN(:,3)))^2);

    for i = 1:nhypos
    AP = [0 0 0] + X(i,:);
    AB = [0 0 0] - [V(1,2) V(2,2) V(3,2)];
    NN(i,:) = [0 0 0] + dot(AP,AB) / dot(AB,AB) * AB;

    end
    W = sqrt((max(NN(:,1))-min(NN(:,1)))^2 + (max(NN(:,2))-min(NN(:,2)))^2 +(max(NN(:,3))-min(NN(:,3)))^2);

end

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








% ex_dist = dot([xb yb zb],AB) / dot(AB,AB) * AB;
% ll= sqrt((ex_dist(1)-min(NN(:,1)))^2 + (ex_dist(2)-min(NN(:,2)))^2 +(ex_dist(3)-min(NN(:,3)))^2);
% lll= sqrt((ex_dist(1)-max(NN(:,1)))^2 + (ex_dist(2)-max(NN(:,2)))^2 +(ex_dist(3)-max(NN(:,3)))^2);
% L = 2*max(ll,lll);

% ex_dist = dot([xb yb zb],AB) / dot(AB,AB) * AB;
% ww= sqrt((ex_dist(1)-min(NN(:,1)))^2 + (ex_dist(2)-min(NN(:,2)))^2 +(ex_dist(3)-min(NN(:,3)))^2);
% www= sqrt((ex_dist(1)-max(NN(:,1)))^2 + (ex_dist(2)-max(NN(:,2)))^2 +(ex_dist(3)-max(NN(:,3)))^2);
% W = 2*max(ww,www);