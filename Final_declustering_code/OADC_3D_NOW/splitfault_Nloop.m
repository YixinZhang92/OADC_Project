i = 0; ii = 0;

while (ii  < 100) && (i< N_loop) 
    ii = ii + 1;

    % Put back the temporary arrays into the original arrays
    % load up arrays with good cluster parameters
    xb = xb_tmp; yb = yb_tmp; zb = zb_tmp; % Barycenters
    xv=xv_tmp; yv=yv_tmp; zv=zv_tmp; % fault plane vertices
    xt=xt_tmp; yt=yt_tmp; zt=zt_tmp; % hypocenter location in a cluster
    vec_plane=vec_plane_tmp; % eigenvector that describes each plane
    Nt=Nt_tmp; % number of events in each trial cluster
    lambda3=lambda3_tmp; % minimum eigenvalue
    L=L_tmp; W=W_tmp; Strike=Strike_tmp; Dip=Dip_tmp; % fault plane parameters

    splitfault(Kfaults)

    % Store each parameter to know which to advance to the next
    % iteration. load up arrays with good cluster parameters
    xb_tmp_ii=xb; yb_tmp_ii=yb; zb_tmp_ii=zb; % Barycenters
    xv_tmp_ii=xv; yv_tmp_ii=yv; zv_tmp_ii=zv; % fault plane vertices
   % xt_tmp_i(:,:,i)=xt; yt_tmp_i(:,:,i)=yt;
   % zt_tmp_i(:,:,i)=zt; hypocenter location in a cluster
    vec_plane_tmp_ii=vec_plane; % eigenvector that describes each plane
    Nt_tmp_ii=Nt; % number of events in each trial cluster
    lambda3_tmp_ii=lambda3; % minimum eigenvalue
    L_tmp_ii=L; W_tmp_ii=W; Strike_tmp_ii=Strike;
    Dip_tmp_ii=Dip; % fault plane parameters

    % increase the fault number
    Kfaults=Kfaults+1;

    if use_glo_var == 1

        JFINALL=faultcluster(con_tol,Kfaults);
    else

        JFINALL=Copy_of_faultcluster(con_tol,Kfaults);
    end

    % Reduce the number of fault because faultcluster.m have
    % increased it by 1.
    Kfaults=Kfaults-1;          

    minDip = min(Dip(1:(Kfaults+1)));
    minNt = min(Nt(1:(Kfaults+1)));

    minDip;
    minNt;

    if (minDip >= dip_threshold) && (i < N_loop) && (minNt >= 4) 
        % minDip >= 0 is normal as if no minDip constraint and no while loop
        i = i+1;

        JFINAL(i) = JFINALL;

        % Store each parameter to know which to advance to the next
        % iteration. load up arrays with good cluster parameters
        xb_tmp_i(:,:,i)=xb_tmp_ii; yb_tmp_i(:,:,i)=yb_tmp_ii; 
        zb_tmp_i(:,:,i)=zb_tmp_ii; % Barycenters
        xv_tmp_i(:,:,i)=xv_tmp_ii; yv_tmp_i(:,:,i)=yv_tmp_ii; 
        zv_tmp_i(:,:,i)=zv_tmp_ii; % fault plane vertices
        % xt_tmp_i(:,:,i)=xt; yt_tmp_i(:,:,i)=yt;
        % zt_tmp_i(:,:,i)=zt; hypocenter location in a cluster
        vec_plane_tmp_i(:,:,i)=vec_plane_tmp_ii; % eigenvector that describes each plane
        Nt_tmp_i(:,:,i)=Nt_tmp_ii; % number of events in each trial cluster
        lambda3_tmp_i(:,:,i)=lambda3_tmp_ii; % minimum eigenvalue
        L_tmp_i(:,:,i)=L_tmp_ii; W_tmp_i(:,:,i)=W_tmp_ii; Strike_tmp_i(:,:,i)=Strike_tmp_ii;
        Dip_tmp_i(:,:,i)=Dip_tmp_ii; % fault plane parameters

        perc = (i/N_loop)*100;
        textprogressbar(perc);
        pause(0.1);

    end                    
end

            
            
            
            
            
            
            
            
%                 % Store each parameter to know which to advance to the
%                 next % iteration. % load up arrays with good cluster
%                 parameters xb_tmp_i(:,:,i)=xb; yb_tmp_i(:,:,i)=yb;
%                 zb_tmp_i(:,:,i)=zb; % Barycenters xv_tmp_i(:,:,i)=xv;
%                 yv_tmp_i(:,:,i)=yv; zv_tmp_i(:,:,i)=zv; % fault plane
%                 vertices xt_tmp_i(:,:,i)=xt; yt_tmp_i(:,:,i)=yt;
%                 zt_tmp_i(:,:,i)=zt;
%                   hypocenter location in a cluster
%                 vec_plane_tmp_i(:,:,i)=vec_plane; % eigenvector that
%                 describes each plane Nt_tmp_i(:,:,i)=Nt; % number of
%                 events in each trial cluster
%                 lambda3_tmp_i(:,:,i)=lambda3; % minimum eigenvalue
%                 L_tmp_i(:,:,i)=L; W_tmp_i(:,:,i)=W;
%                 Strike_tmp_i(:,:,i)=Strike;

 
 
 
 
%   min(Dip_tmp_i(:,1:(Kfaults+1),i))
%    min(Dip_tmp_i(:,1:(Kfaults+1),i))

