function classifying_clusters_from_OADC_2D()

global vec_plane orig_xs orig_ys orig_zs xs ys
global xt yt Nt xb yb lambda3
global L xv yv fscale
global Strike FM_file dist2FM_threshold dip_threshold N_thresh
global xb_tmp_i yb_tmp_i
global xv_tmp_i yv_tmp_i
global xt_tmp_i yt_tmp_i
global vec_plane_tmp_i 
global Nt_tmp_i lambda3_tmp_i
global L_tmp_i Strike_tmp_i err_av kmin kmax N_loop simul_tag infile
global index use_glo_var con_tol Kfaults database database_lambda_only
global line_dens_incr 

num_of_clus = Kfaults;
nclus_now = 0;

for iii=1:num_of_clus

    % Get info of the hypos in the cluster
    xst = xt(iii,1:Nt(iii)); 
    yst = yt(iii,1:Nt(iii));

    xv_clus = xv(iii,:);
    yv_clus = yv(iii,:);

    Lst = L(iii);
    vec_planest = vec_plane(iii,:);

    % determining the variance of each cluster.
    %J_clus=calc_J_per_clus(xst,yst,xv_clus,yv_clus,Lst,vec_planest);
    J_clus = lambda3(iii);             

    if Nt(iii) >= N_thresh && J_clus <= err_av 

        % Simple line density test                
        % determine the density of points along the best-fit line
        xyst = [xst' yst'];

        % compute the covariance matrix for this cluster
        Cxy=cov(xyst,0);

        % Seun checks if Cxy contains NaN.
        NrNaN = sum(isnan(Cxy(:)));
        if NrNaN > 0
            %continue
        end

        % compute the eigenvalues and eigenvectors for this cluster
        [V,D]=eig(Cxy);

        % plot density 
        projX=[V(:,2) V(:,1)]\xyst';
        pxs = projX(1,:);
        pys = projX(2,:);
%
%                 figure
%                 plot(xst,yst,'ro'); hold on;
%                 plot(-pxs,zeros(1,length(pxs)),'bo'); shg

        % line density will be equal to 0 if the cluster pass the line density test.
        % and 1 if it doesn't.
        %figure;
        %h = histogram(pxs, min(pxs):max(pxs)+1);
        %hist_values = h.Values(1:end-1);
        %line_density = max(hist_values < N_thresh);

        hist_values = histcounts(pxs, min(pxs):line_dens_incr:max(pxs)+line_dens_incr);
        line_density = max(hist_values(1:end-1) < 1); %
        %line_density=0;

        if line_density == 0
            nclus_now = nclus_now + 1;

            for ehypo=1:Nt(iii) % for each hypo in the cluster
                % get the index of each hypocenter in the original catalog
                index_hypo = find(xs==xt(iii,ehypo) & ys==yt(iii,ehypo));

                % assign J_clus to each hypocenter (using lambda2 only)
                if (J_clus < database_lambda_only(index_hypo,4))
                    database_lambda_only(index_hypo,4) = J_clus;
                    database_lambda_only(index_hypo,5) = Nt(iii);
                end
                
                % assign N and J_clus to each hypocenter (using both lambda2 and Neqs)
                if (Nt(iii) > database(index_hypo,5))
                    database(index_hypo,4) = J_clus;
                    database(index_hypo,5) = Nt(iii);
                end
                
            end
        end    
    end
end
% end of hypocenter classification
