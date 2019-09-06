function fit_planes_and_plot_clusters_based_on_lambda2_and_Neqs()

global  N_thresh simul_tag database
global L_ln W_ln Strike_ln Dip_ln xv_ln yv_ln zv_ln vec_plane_ln lambda3_ln

colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'c.', 'm.'};

unq_nclus = unique([database(:,4) database(:,5)],'rows');
[tot_num_clus,~]= size(unq_nclus); kk = 0;

for i=1:tot_num_clus
    
    xst = database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),1); 
    yst = database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),2);
    zst = database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),3);
    
    if (length(xst)>=N_thresh) && (unq_nclus(i,1)~=100)
        % compute the covariance matrix for this cluster
        Cxy=cov([xst yst zst],0);

        % Seun checks if Cxy contains NaN.
        NrNaN = sum(isnan(Cxy(:)));
        if NrNaN > 0
            continue
        end

        kk = kk + 1;
        xbt(kk) = mean(xst); ybt(kk) = mean(yst); zbt(kk) = mean(zst);

        % compute the eigenvalues and eigenvectors for this cluster
        [V,D]=eig(Cxy);

        % calculate fault plane parameters from the eigen results
        % and calculate the vertices of the fault plane
        XX = [xst yst zst];

        [L_ln(kk),W_ln(kk),Strike_ln(kk),Dip_ln(kk),xv_ln(kk,:),yv_ln(kk,:),zv_ln(kk,:)] = ...
            fltplane(XX,V,D,xbt(kk),ybt(kk),zbt(kk));

        % save the plane unit normal vector and eigenvalue
        vec_plane_ln(kk,1:3)=V(1:3,1);
        lambda3_ln(kk)=sqrt(D(1,1));
    end
end

% Ploting figures
Fig1 = figure('Name','Clustered hypocenters with fault model','Position', get(0, 'Screensize'));

for i=1:tot_num_clus  
    if unq_nclus(i,1) ~= 100
        plot3(database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),1),...
            database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),2),...
            database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),3),...
            colors{mod(i,7)+1},'MarkerSize',12); hold on
    else
        plot3(database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),1),...
            database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),2),...
            database(database(:,4)==unq_nclus(i,1) & database(:,5)==unq_nclus(i,2),3),'k.'); hold on            
    end
end

if kk>0
    for m=1:kk
        fill3(xv_ln(m,1:4),yv_ln(m,1:4),zv_ln(m,1:4),...
            'w','FaceAlpha',0.2,'FaceColor',[0.5 0.5 0.5]);
    end
end

grid on; axis equal;
title('Fault Model (Constraint: Lambda2 and Neqs)');
set(gca, 'fontsize', 18); shg

% Printing figure to file
fig_filename = [simul_tag '.faultmodel.lambda2.and.Neqs.png'];
F1    = getframe(Fig1);
imwrite(F1.cdata, fig_filename, 'png')
savefig(Fig1,[fig_filename(1:end-4) '.fig']);



% Fig1 = figure('Name','Clustered hypocenters','Position', get(0, 'Screensize'));
% 
% for i=1:length(unq_nclus)  
%     if unq_nclus(i) ~= 100
%         plot3(database(database(:,4)==unq_nclus(i),1),database(database(:,4)==unq_nclus(i),2),...
%             database(database(:,4)==unq_nclus(i),3),colors{mod(i,7)+1},'MarkerSize',12); hold on
%     else
%         plot3(database(database(:,4)==unq_nclus(i),1),database(database(:,4)==unq_nclus(i),2),...
%             database(database(:,4)==unq_nclus(i),3),'k.'); hold on            
%     end
% end
% 
% set(gca, 'fontsize', 18);
% 
% % Printing figure to file
% fig_filename = [simul_tag '.final.model.png'];
% F1    = getframe(Fig1);
% imwrite(F1.cdata, fig_filename, 'png')
% savefig(Fig1,[fig_filename(1:end-4) '.fig'])



% function fit_planes_and_plot_clusters_based_on_lambda2_and_Neqs()
% 
% global  N_thresh simul_tag database
% global L_ln W_ln Strike_ln Dip_ln xv_ln yv_ln zv_ln vec_plane_ln lambda3_ln
% 
% colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'c.', 'm.'};
% 
% unq_nclus = unique(database(:,4));
% tot_num_clus= length(unq_nclus); kk = 0;
% 
% for i=1:tot_num_clus
%     
%     xst = database(database(:,4)==unq_nclus(i),1); 
%     yst = database(database(:,4)==unq_nclus(i),2);
%     zst = database(database(:,4)==unq_nclus(i),3);
%     
%     if (length(xst)>=N_thresh) && (unq_nclus(i)~=100)
%         % compute the covariance matrix for this cluster
%         Cxy=cov([xst yst zst],0);
% 
%         % Seun checks if Cxy contains NaN.
%         NrNaN = sum(isnan(Cxy(:)));
%         if NrNaN > 0
%             continue
%         end
% 
%         kk = kk + 1;
%         xbt(kk) = mean(xst); ybt(kk) = mean(yst); zbt(kk) = mean(zst);
% 
%         % compute the eigenvalues and eigenvectors for this cluster
%         [V,D]=eig(Cxy);
% 
%         % calculate fault plane parameters from the eigen results
%         % and calculate the vertices of the fault plane
%         XX = [xst yst zst];
% 
%         [L_ln(kk),W_ln(kk),Strike_ln(kk),Dip_ln(kk),xv_ln(kk,:),yv_ln(kk,:),zv_ln(kk,:)] = ...
%             fltplane(XX,V,D,xbt(kk),ybt(kk),zbt(kk));
% 
%         % save the plane unit normal vector and eigenvalue
%         vec_plane_ln(kk,1:3)=V(1:3,1);
%         lambda3_ln(kk)=sqrt(D(1,1));
%     end
% end
% 
% % Ploting figures
% Fig1 = figure('Name','Clustered hypocenters with fault model','Position', get(0, 'Screensize'));
% 
% for i=1:length(unq_nclus)  
%     if unq_nclus(i) ~= 100
%         plot3(database(database(:,4)==unq_nclus(i),1),database(database(:,4)==unq_nclus(i),2),...
%             database(database(:,4)==unq_nclus(i),3),colors{mod(i,7)+1},'MarkerSize',12); hold on
%     else
%         plot3(database(database(:,4)==unq_nclus(i),1),database(database(:,4)==unq_nclus(i),2),...
%             database(database(:,4)==unq_nclus(i),3),'k.'); hold on            
%     end
% end
% 
% if kk>0
%     for m=1:kk
%         fill3(xv_ln(m,1:4),yv_ln(m,1:4),zv_ln(m,1:4),...
%             'w','FaceAlpha',0.2,'FaceColor',[0.5 0.5 0.5]);
%     end
% end
% 
% grid on; axis equal;
% title('Fault Model (Constraint: Lambda2 and Neqs)');
% set(gca, 'fontsize', 18); shg
% 
% % Printing figure to file
% fig_filename = [simul_tag '.faultmodel.lambda2.and.Neqs.png'];
% F1    = getframe(Fig1);
% imwrite(F1.cdata, fig_filename, 'png')
% savefig(Fig1,[fig_filename(1:end-4) '.fig']);
% 
% 
% 
% % Fig1 = figure('Name','Clustered hypocenters','Position', get(0, 'Screensize'));
% % 
% % for i=1:length(unq_nclus)  
% %     if unq_nclus(i) ~= 100
% %         plot3(database(database(:,4)==unq_nclus(i),1),database(database(:,4)==unq_nclus(i),2),...
% %             database(database(:,4)==unq_nclus(i),3),colors{mod(i,7)+1},'MarkerSize',12); hold on
% %     else
% %         plot3(database(database(:,4)==unq_nclus(i),1),database(database(:,4)==unq_nclus(i),2),...
% %             database(database(:,4)==unq_nclus(i),3),'k.'); hold on            
% %     end
% % end
% % 
% % set(gca, 'fontsize', 18);
% % 
% % % Printing figure to file
% % fig_filename = [simul_tag '.final.model.png'];
% % F1    = getframe(Fig1);
% % imwrite(F1.cdata, fig_filename, 'png')
% % savefig(Fig1,[fig_filename(1:end-4) '.fig'])
% 
