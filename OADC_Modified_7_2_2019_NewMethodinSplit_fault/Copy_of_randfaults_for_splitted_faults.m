function Copy_of_randfaults_for_splitted_faults(FAULT_FLAG)

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
%FAULT_FLAG

% vertices of the thick fault
xv1 = xv(FAULT_FLAG,1); xv2 = xv(FAULT_FLAG,2); 
xv3 = xv(FAULT_FLAG,3); xv4 = xv(FAULT_FLAG,4); 

yv1 = yv(FAULT_FLAG,1); yv2 = yv(FAULT_FLAG,2); 
yv3 = yv(FAULT_FLAG,3); yv4 = yv(FAULT_FLAG,4); 

zv1 = zv(FAULT_FLAG,1); zv2 = zv(FAULT_FLAG,2); 
zv3 = zv(FAULT_FLAG,3); zv4 = zv(FAULT_FLAG,4); 
    
L2=L(FAULT_FLAG)./2.0;
W2=W(FAULT_FLAG);

for k=1:2
    
    % Seun did not make the L to be greater than W, since the vertices
    % will be determined from the vertices of the thick fault.
    L(k)=L2;
    W(k)=W2;
        
%         
%     %  make L2 >= W2
%     if L2 <= W2
%         L(k)=W2;
%         W(k)=L2;
%     else
%         L(k)=L2;
%         W(k)=W2;
%     end
    
    % Seun changed the algorithm to preserve the orientation of the original
    % fault. So, the two new faults have the same orientation but at random
    % positions.

    vec_plane(k,1:3)=vec_plane(FAULT_FLAG,:);
    
    % now randomly pick a hypocenter location to be the center of the fault
    nb=randperm(Nt(FAULT_FLAG));
    xb(k)=xt(FAULT_FLAG,nb(1));
    yb(k)=yt(FAULT_FLAG,nb(1));
    zb(k)=zt(FAULT_FLAG,nb(1));
    

    
    % compute vertice locations of the two subfaults
    if k==1
        % first subfault
        % Middle point
        xr1 = (3*xv1 + xv2 + xv3 +3*xv4)/8;
        yr1 = (3*yv1 + yv2 + yv3 +3*yv4)/8;
        zr1 = (3*zv1 + zv2 + zv3 +3*zv4)/8;
        
        % Uncomment to apply no displacement to subfault 1
%         xb(k) = xr1;
%         yb(k) = yr1;
%         zb(k) = zr1;

    
        % correction factor to move the new coordinate to the random
        % position as its centroid
        phi_xb = xb(k) - xr1;
        phi_yb = yb(k) - yr1;
        phi_zb = zb(k) - zr1;
        
       
        xvt(k,1)= xv1 + phi_xb;
        yvt(k,1)= yv1 + phi_yb;
        zvt(k,1)= zv1 + phi_zb;

        xvt(k,2)= (xv1 + xv2)/2 + phi_xb;
        yvt(k,2)= (yv1 + yv2)/2 + phi_yb;
        zvt(k,2)= (zv1 + zv2)/2 + phi_zb;

        xvt(k,3)= (xv3 + xv4)/2 + phi_xb;
        yvt(k,3)= (yv3 + yv4)/2 + phi_yb;
        zvt(k,3)= (zv3 + zv4)/2 + phi_zb;

        xvt(k,4)= xv4 + phi_xb;
        yvt(k,4)= yv4 + phi_yb;
        zvt(k,4)= zv4 + phi_zb;
        
        
    else
        % second subfault
        % Middle point
        xr2 = (xv1 + 3*xv2 + 3*xv3 +xv4)/8;
        yr2 = (yv1 + 3*yv2 + 3*yv3 +yv4)/8;
        zr2 = (zv1 + 3*zv2 + 3*zv3 +zv4)/8;
        
        % Uncomment to apply no displacement to subfault 2
%         xb(k) = xr2;
%         yb(k) = yr2;
%         zb(k) = zr2;

        
        % correction factor to move the new coordinate to the random
        % position as its centroid
        phi_xb = xb(k)- xr2;
        phi_yb = yb(k)- yr2;
        phi_zb = zb(k) - zr2;
        
        
        xvt(k,1)= (xv1 + xv2)/2 + phi_xb;
        yvt(k,1)= (yv1 + yv2)/2 + phi_yb;
        zvt(k,1)= (zv1 + zv2)/2 + phi_zb;

        xvt(k,2)= xv2 + phi_xb;
        yvt(k,2)= yv2 + phi_yb;
        zvt(k,2)= zv2 + phi_zb;

        xvt(k,3)= xv3 + phi_xb;
        yvt(k,3)= yv3 + phi_yb;
        zvt(k,3)= zv3 + phi_zb;

        xvt(k,4)= (xv3 + xv4)/2 + phi_xb;
        yvt(k,4)= (yv3 + yv4)/2 + phi_yb;
        zvt(k,4)= (zv3 + zv4)/2 + phi_zb;
        
        
    end
end

xv = xvt;
yv = yvt;
zv = zvt;


% Uncomment to display thich fault and its subfaults
%
% figure
% plot(xv1,yv1,'ro'); hold on;
% plot(xv2,yv2,'bo'); hold on;
% plot(xv3,yv3,'go'); hold on;
% plot(xv4,yv4,'ko'); hold on;
% plot([xv1,xv2,xv3,xv4,xv1],[yv1,yv2,yv3,yv4,yv1],'g','LineWidth',3.5); hold on;
% 
% %subfault1
% k = 1;
% plot(xvt(k,1),yvt(k,1),'ro'); hold on;
% plot(xvt(k,2),yvt(k,2),'bo'); hold on;
% plot(xvt(k,3),yvt(k,3),'go'); hold on;
% plot(xvt(k,4),yvt(k,4),'ko'); hold on;
% plot([xvt(k,1),xvt(k,2),xvt(k,3),xvt(k,4),xvt(k,1)],[yvt(k,1),yvt(k,2),yvt(k,3),yvt(k,4),yvt(k,1)], 'r','LineWidth',2); hold on;
% 
% 
% %subfault2
% k = 2;
% plot(xvt(k,1),yvt(k,1),'ro'); hold on;
% plot(xvt(k,2),yvt(k,2),'bo'); hold on;
% plot(xvt(k,3),yvt(k,3),'go'); hold on;
% plot(xvt(k,4),yvt(k,4),'ko'); hold on;
% plot([xvt(k,1),xvt(k,2),xvt(k,3),xvt(k,4),xvt(k,1)],[yvt(k,1),yvt(k,2),yvt(k,3),yvt(k,4),yvt(k,1)], 'b','LineWidth',2); hold on;

return;
   