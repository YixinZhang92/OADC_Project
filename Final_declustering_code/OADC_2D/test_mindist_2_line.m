clear all; close all; clc;
% the minimum distance from a 2D point to a line in 2D

xs = [3 7 11 8 5 0];
ys = [1 5 18 23 18 5];
xv = [3 8];
yv = [4 18.1421];

plot(xv,yv,'r'); xlim([0 15]); ylim([0 25]); hold on; 
plot(xs,ys,'bo'); grid MINOR

for k=1:6
    % Choose hypocenter vector
    u=[xs(k) ys(k);xs(k) ys(k)];
    u1=[xs(k) ys(k)];
    
    % Construct vertices vector
    v=[ xv' yv'];
    
    % Construct distance vector from each vertex
    d=u-v;
    d1=[d(1,1) d(1,2)];
    d2=[d(2,1) d(2,2)];
    d1n=norm(d1);
    d2n=norm(d2);
    L1=sqrt((xv(2)-xv(1))^2+(yv(2)-yv(1))^2);
    
    % Find all angles and perpendicular distances to edge
    % 12 Face
    %gam12=acos(d12./(d1n.*d2n))
    [alph12,beta12,gam12]=trisol(d2n,d1n,L1,'r');
    a12=d1n.*d2n.*sin(gam12)./L1;
    
    % compute the covariance matrix for this cluster
    Cxy=cov( [xv' yv'],0);
    
    % compute the eigenvalues and eigenvectors for this cluster
    [V,D]=eig(Cxy);
    
    % save the plane unit normal vector and eigenvalue
    vec_plane(1:2)=V(1:2,1);
    xb = mean(xv);
    yb = mean(yv);
    
    % Find the perpendicular distance to the infinite plane
    vec(1:2)=vec_plane(1:2);
    u2=u1-[xb yb];
    dperp=abs(dot(u2,vec));
    
    % Find the position field for the hypocenter relative to the plane and find
    % the minimum distance.
    p2=pi/2.0;
    dmin(k)=dperp;
    if (alph12 > p2); dmin(k)=d1n; end
    if (beta12 > p2); dmin(k)=d2n; end
    
end
dmin