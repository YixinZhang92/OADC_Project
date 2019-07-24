function Copy_of_randfaults(n0,FAULT_FLAG)

%  Construct n0 randomn faults
%  just need the vertice locations

%  FAULT_FLAG = 0,  Use all hypocenters.  Initialization of the Process
%  FAULT_FLAG = kthick (from splitfault.m),  Split the thickest fault into
%               two random planes of length L/2

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt yt zt Nt xb yb zb lambda3

%fprintf('From randfaults\n');

for k=1:n0;
    
    % generate a random fault using a random square 3x3 matrix.  The
    % smallest eigenvalue will be the fault normal.
    A=rand(3);%eye(3);%
    B=A'*A;
    [V,D]=eig(B);
    
    % choose 2 random numbers between 0 and 1
    L1=rand(1,2);
    % choose a permutation of 1-10
    a=randperm(10);
    
    if FAULT_FLAG == 0;
        L2=L1(1).*fscale;
        W2=L1(2).*fscale;
    else
        L2=L(FAULT_FLAG)./2.0;
        W2=W(FAULT_FLAG);
    end
    
    %  make L2 >= W2
    if L2 <= W2;
        L(k)=W2;
        W(k)=L2;
    else;
        L(k)=L2;
        W(k)=W2;
    end;
    
    %  find the plane unit vector
    vec_plane(k,1:3)=V(1:3,1);
    
    % now randomly pick a hypocenter location to be the center of the fault
    if FAULT_FLAG == 0;
        nb=randperm(N);
        xb(k)=xs(nb(1));%mean(xs);%
        yb(k)=ys(nb(1));%mean(ys);%
        zb(k)=zs(nb(1));%mean(zs);%
    else
        nb=randperm(Nt(FAULT_FLAG));
        xb(k)=xt(FAULT_FLAG,nb(1));
        yb(k)=yt(FAULT_FLAG,nb(1));
        zb(k)=zt(FAULT_FLAG,nb(1));
    end

    
    L2=L(k)./2.0;
    W2=W(k)./2.0;

    % compute vertice locations
    xv(k,1)=W2.*V(1,2)+L2.*V(1,3) + xb(k);
    yv(k,1)=W2.*V(2,2)+L2.*V(2,3) + yb(k);
    zv(k,1)=W2.*V(3,2)+L2.*V(3,3) + zb(k);

    xv(k,2)=W2.*V(1,2)-L2.*V(1,3) + xb(k);
    yv(k,2)=W2.*V(2,2)-L2.*V(2,3) + yb(k);
    zv(k,2)=W2.*V(3,2)-L2.*V(3,3) + zb(k);

    xv(k,3)=-W2.*V(1,2)-L2.*V(1,3) + xb(k);
    yv(k,3)=-W2.*V(2,2)-L2.*V(2,3) + yb(k);
    zv(k,3)=-W2.*V(3,2)-L2.*V(3,3) + zb(k);

    xv(k,4)=-W2.*V(1,2)+L2.*V(1,3) + xb(k);
    yv(k,4)=-W2.*V(2,2)+L2.*V(2,3) + yb(k);
    zv(k,4)=-W2.*V(3,2)+L2.*V(3,3) + zb(k);
    
end;

return;
    
    