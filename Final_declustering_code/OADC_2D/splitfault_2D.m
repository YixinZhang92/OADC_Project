function splitfault_2D(Kfaults)
%  splitfault - split the fault with the largest lambda3 eigenvalue into 2.

% n0 = number of fault clusters to determine
% J  = global variance is output

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

% Seun's edits: % The cluster with the greatest fault thickness based on the minimum
% eigenvalue will now be split into two parts. If a file containign focal 
% mechanisms is available, the focal mechanisms will be used to determine 
% the orientation of the two faults instead of using the random-seeded planes. 


global xc yc vec_plane xb_old yb_old xs ys N Nc
global xt yt Nt xb yb lambda3
global L xv yv L_old xv_old yv_old fscale
global xt_old yt_old vec_plane_old lambda3_old
global Strike Strike_old Nt_old

%  find the index of the cluster with the largest lambda3
[lambda3max,kthick]=max(lambda3);

%  Save the fault models for the better fitting faults
kg=0;
for k=1:Kfaults
    
    if k ~= kthick
        kg=kg+1;
        % load up arrays with good cluster parameters
        
        % Barycenters
        xb_old(kg)=xb(k);
        yb_old(kg)=yb(k);        
        
        % fault plane vertices
        xv_old(kg,:)=xv(k,:);
        yv_old(kg,:)=yv(k,:);
        
        % hypocenter location in a cluster
        xt_old(kg,:)=xt(k,:);
        yt_old(kg,:)=yt(k,:);
        
        % eigenvector that describes each plane
        vec_plane_old(kg,1:2)=vec_plane(k,:);

        % number of events in each trial cluster
        Nt_old(kg)=Nt(k);

        % minimum eigenvalue
        lambda3_old(kg)=lambda3(k);

        % fault plane parameters
        L_old(kg)=L(k);
        Strike_old(kg)=Strike(k);
       
    end
end

% % The cluster with the greatest fault thickness based on the minimum
% % eigenvalue will now be split into two random parts. This is the "kthick"
% % cluster
% 
randfaults_2D(2,kthick);

% The cluster with the greatest fault thickness based on the minimum
% eigenvalue will now be split into two parts. If a file containign focal 
% mechanisms is available, the focal mechanisms will be used to determine 
% the orientation of the two faults instead of using the random-seeded planes. 
% This is the "kthick" cluster.

%FM_seeded_planes(2,kthick)

%  Now add these new faults to the other clusters in the "old" storage
for k=1:2 
        kg=kg+1;
        % load up old arrays with new cluster parameters
        
        % Barycenters
        xb_old(kg)=xb(k);
        yb_old(kg)=yb(k);
        
        % fault plane vertices
        xv_old(kg,:)=xv(k,:);
        yv_old(kg,:)=yv(k,:);
        
        % eigenvector that describes each plane
        vec_plane_old(kg,1:2)=vec_plane(k,:);

        % fault plane parameters
        L_old(kg)=L(k);
        
end

%  Load up working arrays with all data

xb=xb_old;
yb=yb_old;

xv=xv_old;
yv=yv_old;

vec_plane=vec_plane_old;
L=L_old;

%  All done.  Start the analysis again

end

