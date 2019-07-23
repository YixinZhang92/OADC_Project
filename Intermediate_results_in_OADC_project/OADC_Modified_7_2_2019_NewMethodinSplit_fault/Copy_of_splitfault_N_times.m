function Copy_of_splitfault_N_times(con_tol, Kfaults)
%  splitfault - split the fault with the largest lambda3 eigenvalue into 2.

% n0 = number of fault clusters to determine
% J  = global variance is output

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt_old yt_old zt_old vec_plane_old lambda3_old
global Strike Dip Strike_old Dip_old Nt_old

%  find the index of the cluster with the largest lambda3
[lambda3max,kthick]=max(lambda3);

%  Save the fault models for the better fitting faults
kg=0;

for k=1:Kfaults
        
    if Kfaults == 1
        kg=kg+1;
        
        % load up arrays with good cluster parameters
        
        % Barycenters
        xb_old(kg)=xb(k);
        yb_old(kg)=yb(k);
        zb_old(kg)=zb(k);
        
        % fault plane vertices
        xv_old(kg,:)=xv(k,:);
        yv_old(kg,:)=yv(k,:);
        zv_old(kg,:)=zv(k,:);
        
        % hypocenter location in a cluster
        xt_old(kg,:)=xt(k,:);
        yt_old(kg,:)=yt(k,:);
        zt_old(kg,:)=zt(k,:);
        
        % eigenvector that describes each plane
        vec_plane_old(kg,1:3)=vec_plane(k,:);

        % number of events in each trial cluster
        Nt_old(kg)=Nt(k);

        % minimum eigenvalue
        lambda3_old(kg)=lambda3(k);

        % fault plane parameters
        L_old(kg)=L(k);
        W_old(kg)=W(k);
        Strike_old(kg)=Strike(k);
        Dip_old(kg)=Dip(k);
       
    else
        
        if k ~= kthick
            kg=kg+1;

            % load up arrays with good cluster parameters

            % Barycenters
            xb_old(kg)=xb(k);
            yb_old(kg)=yb(k);
            zb_old(kg)=zb(k);

            % fault plane vertices
            xv_old(kg,:)=xv(k,:);
            yv_old(kg,:)=yv(k,:);
            zv_old(kg,:)=zv(k,:);

            % hypocenter location in a cluster
            xt_old(kg,:)=xt(k,:);
            yt_old(kg,:)=yt(k,:);
            zt_old(kg,:)=zt(k,:);

            % eigenvector that describes each plane
            vec_plane_old(kg,1:3)=vec_plane(k,:);

            % number of events in each trial cluster
            Nt_old(kg)=Nt(k);

            % minimum eigenvalue
            lambda3_old(kg)=lambda3(k);

            % fault plane parameters
            L_old(kg)=L(k);
            W_old(kg)=W(k);
            Strike_old(kg)=Strike(k);
            Dip_old(kg)=Dip(k);
        end
    end
end

% The cluster with the greatest fault thickness based on the minimum
% eigenvalue will now be split into two random parts. This is the "kthick"
% cluster

%% Temporary store the parameters somewhereelse
% % load up arrays with good cluster parameters
% 
% % Barycenters
% xb_old_tmp = xb_old;
% yb_old_tmp = yb_old;
% zb_old_tmp = zb_old;
% 
% % fault plane vertices
% xv_old_tmp = xv_old;
% yv_old_tmp = yv_old;
% zv_old_tmp = zv_old;
% 
% % hypocenter location in a cluster
% xt_old_tmp = xt_old;
% yt_old_tmp = yt_old;
% zt_old_tmp = zt_old;
% 
% % eigenvector that describes each plane
% vec_plane_old_tmp = vec_plane_old;
% 
% % number of events in each trial cluster
% Nt_old_tmp = Nt_old;
% 
% % minimum eigenvalue
% lambda3_old_tmp = lambda3_old;
% 
% % fault plane parameters
% L_old_tmp = L_old;
% W_old_tmp = W_old;
% Strike_old_tmp = Strike_old;
% Dip_old_tmp = Dip_old;

%% Loop over the split N times, and find the configuration that has the lowest lambda3

n_loop = 5;
kg
for i = 1:n_loop % ********* loop N times *********
        
        % store kg somewhere to avoid overwritten
        kg_split = kg;

        % Split the thickest fault into two faults
        Copy_of_randfaults(2,kthick);
        % Seun changed the algorithm to preserve the orientation of the original
        % fault. So, the two new faults have the same orientation but at random
        % positions.
        %Copy_of_randfaults_for_splitted_faults(kthick)

        %  Now add these new faults to the other clusters in the "old" storage

        for k=1:2
                if (Kfaults == 1) && k ==1
                    kg_split =1;
                elseif (Kfaults == 1) && k ==2
                    kg_split = 2;
                else
                    kg_split=kg_split+1;
                end
                
                kg_split
                % load up old arrays with new cluster parameters

                % Barycenters
                xb_old(kg_split)=xb(k);
                yb_old(kg_split)=yb(k);
                zb_old(kg_split)=zb(k);

                % fault plane vertices
                xv_old(kg_split,:)=xv(k,:);
                yv_old(kg_split,:)=yv(k,:);
                zv_old(kg_split,:)=zv(k,:);

                % eigenvector that describes each plane
                vec_plane_old(kg_split,1:3)=vec_plane(k,:);

                % fault plane parameters
                L_old(kg_split)=L(k);
                W_old(kg_split)=W(k);
                
                
                % temporary store the parameters in this loop 
                % Barycenters
                
                xb_tmp(i,k)=xb(k);
                yb_tmp(i,k)=yb(k);
                zb_tmp(i,k)=zb(k);

                % fault plane vertices
                xv_tmp(i,k,:)=xv(k,:);
                yv_tmp(i,k,:)=yv(k,:);
                zv_tmp(i,k,:)=zv(k,:);

                % eigenvector that describes each plane
                vec_plane_tmp(i,k,:)=vec_plane(k,:);

                % fault plane parameters
                L_tmp(i,k)=L(k);
                W_tmp(i,k)=W(k);

        end

        
        
        
        %  Load up working arrays with all data

        xb=xb_old;
        yb=yb_old;
        zb=zb_old;

        xv=xv_old;
        yv=yv_old;
        zv=zv_old;

        vec_plane=vec_plane_old;
        L=L_old;
        W=W_old;
    
        JFINAL(i)=Copy_of_faultcluster_using_lamda3(con_tol,Kfaults);


end


JFINAL
xb_tmp
%  All done.  Start the analysis again

end

