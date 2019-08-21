function [xbt,ybt,zbt,xvt,yvt,zvt,vec_plane_t,Ltn,Wtn] = randfaults_local(xst,yst,zst,Lt,Wt,n0)



n0=2;

for k=1:n0
    
    % generate a random fault using a random square 3x3 matrix.  The
    % smallest eigenvalue will be the fault normal.
    A=rand(3);
    B=A'*A;
    [V,D]=eig(B);
    
    % choose 2 random numbers between 0 and 1
    L1=rand(1,2);
    % choose a permutation of 1-10
    a=randperm(10);
    
    L2=Lt./2.0;
    W2=Wt;
    
    %  make L2 >= W2
    if L2 <= W2
        Ltn(k)=W2;
        Wtn(k)=L2;
    else
        Ltn(k)=L2;
        Wtn(k)=W2;
    end
    
    %  find the plane unit vector
        % Seun changed the algorithm to preserve the orientation of the original
    % fault. So, the two new faults have the same orientation but at random
    % positions.
   
    vec_plane_t(k,1:3)=V(1:3,1);        

    % now randomly pick a hypocenter location to be the center of the fault
    %FAULT_FLAG
    nb=randperm(length(xst)); 
    xbt(k)=xst(nb(1));
    ybt(k)=yst(nb(1));
    zbt(k)=zst(nb(1));
    
    L2=Ltn(k)./2.0;
    W2=Wtn(k)./2.0;

    % compute vertice locations
    xvt(k,1)=W2.*V(1,2)+L2.*V(1,3) + xb(k);
    yvt(k,1)=W2.*V(2,2)+L2.*V(2,3) + yb(k);
    zvt(k,1)=W2.*V(3,2)+L2.*V(3,3) + zb(k);

    xvt(k,2)=W2.*V(1,2)-L2.*V(1,3) + xb(k);
    yvt(k,2)=W2.*V(2,2)-L2.*V(2,3) + yb(k);
    zvt(k,2)=W2.*V(3,2)-L2.*V(3,3) + zb(k);

    xvt(k,3)=-W2.*V(1,2)-L2.*V(1,3) + xb(k);
    yvt(k,3)=-W2.*V(2,2)-L2.*V(2,3) + yb(k);
    zvt(k,3)=-W2.*V(3,2)-L2.*V(3,3) + zb(k);

    xvt(k,4)=-W2.*V(1,2)+L2.*V(1,3) + xb(k);
    yvt(k,4)=-W2.*V(2,2)+L2.*V(2,3) + yb(k);
    zvt(k,4)=-W2.*V(3,2)+L2.*V(3,3) + zb(k);
    
end


JFINAL=faultcluster_local(xvt,yvt,zvt,xbt,ybt,zbt,vec_plane_t,Ltn,Wtn,con_tol,Kfaults);
    

return;
    
    
