function randfaults_2D(n0,FAULT_FLAG)

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

global xc yc vec_plane xb_old yb_old xs ys N Nc
global L xv yv L_old xv_old yv_old fscale
global xt yt Nt xb yb lambda3

%fprintf('From randfaults\n');

for k=1:n0
    
    % generate a random fault using a random square 3x3 matrix.  The
    % smallest eigenvalue will be the fault normal.
    A=rand(2);
    B=A'*A;
    [V,D]=eig(B);
    
    % choose 2 random numbers between 0 and 1
    L1=rand(1,2);
    % choose a permutation of 1-10
    a=randperm(10);
    
    if FAULT_FLAG == 0
        L2=L1(1).*fscale;
        %W2=L1(2).*fscale;
    else
        L2=L(FAULT_FLAG)./2.0;
        %W2=W(FAULT_FLAG);
    end
    
%     %  make L2 >= W2
%     if L2 <= W2
%         L(k)=W2;
%         W(k)=L2;
%     else
         L(k)=L2;
%         W(k)=W2;
%     end
    
    %  find the plane unit vector
        % Seun changed the algorithm to preserve the orientation of the original
    % fault. So, the two new faults have the same orientation but at random
    % positions.
    
%     if FAULT_FLAG == 0
        %  find the plane unit vector
        vec_plane(k,1:2)=V(1:2,1);        
%     else
%         vec_plane(k,1:3)=vec_plane(FAULT_FLAG,:);
%         V(1:3,1)=vec_plane(FAULT_FLAG,:); 
%     end


    % now randomly pick a hypocenter location to be the center of the fault
    if FAULT_FLAG == 0
        nb=randperm(N);
        xb(k)=xs(nb(1));
        yb(k)=ys(nb(1));
        
    else
        %FAULT_FLAG
        Nt(FAULT_FLAG);
        nb=randperm(Nt(FAULT_FLAG)); 
        xb(k)=xt(FAULT_FLAG,nb(1));
        yb(k)=yt(FAULT_FLAG,nb(1));
        
    end

    
    L2=L(k)./2.0;
   
    % compute vertice locations
    xv(k,1)=L2.*V(1,2) + xb(k);
    yv(k,1)=L2.*V(2,2) + yb(k);

    xv(k,2)=-L2.*V(1,2) + xb(k);
    yv(k,2)=-L2.*V(2,2) + yb(k);
    
end

return;
    
    
