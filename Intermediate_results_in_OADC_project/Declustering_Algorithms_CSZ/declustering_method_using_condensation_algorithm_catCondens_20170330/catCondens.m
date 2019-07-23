function con_cat=catCondens(org_cat,cov_mat,DO_PLOT)
% A function for condensation of location distributions for optimal spatial information encoding
% 
% For detailed information check the following publication:
% Y. Kamer, G. Ouillon, D. Sornette and J. Wössner (2014)
% Condensation of earthquake location distributions: Optimal spatial information encoding 
% and application to multifractal % analysis of South Californian seismicity
% http://arxiv.org/abs/1407.1036
% 
% Function syntax
% org_cat : N by D matrix containing coordinates of N points in D dimensions. D can be 2 or 3.
% cov_mat : N by C matrix containing location uncertainty for each event
%   if D=C=2 || D=C=3 : expects standard error (standard deviation) for each dimension
%   if D=2  C=3: expects elements covariance matrix covXX,covYY,covXY
%   if D=3  C=6: expects elements covariance matrix covXX,covYY,covZZ,covXY,covXZ,covYZ
% DO_PLOT: if true generates scatter plots of the original and condensed catalogs. 
%   Colormaps are saturated at 85th percentile of weight/ isotropic variance
%
% Y.Kamer, 
% Zurich 20141209
%%%%%
%
% VERSION HISTORY
% 20170330: Added the randnorm.m function
% 20150721: Added 2D/3D ellipsoid plots using Gautam Vallabha's function
% 20150410: Corrected a bug which caused some the off-diagonal elements in the covariance matrices to be set to zero
% 
%%%%%
% Example : 200 uniformly distributed points
% org_cat=rand(200,2)*20;
% cov_mat=rand(200,2)*3;
% con_cat=catCondens(org_cat,cov_mat,1);

if(nargin==2)
    DO_PLOT=0;
end
[N,D]   = size(org_cat);
C       = size(cov_mat,2);
numMO   = 100; %Number of sample points for the sampling of the source pdfs 

    cov_matTMP = zeros(N,D^2);
    if(C==D)
        cov_matTMP(:,1:D+1:D^2)     = cov_mat.^2;
    elseif(D==2 && C==3)
        cov_matTMP(:,[1 4 3])       = cov_mat;
        cov_matTMP(:,2)             = cov_matTMP(:,3);
    elseif(D==3 && C==6)
        cov_matTMP(:,[1 5 9 4 7 8]) = cov_mat;
        cov_matTMP(:,[2 3 6])       = cov_matTMP(:,[4 7 8]);
    else
        disp('Inconsistent number of columns in the covariance matrix!')
        return;
    end
    cov_mat = cov_matTMP;
    con_cat  = [org_cat ones(size(org_cat,1),1)];
    
    tot_var      = zeros(N,1);
    wght_vec     = con_cat(:,D+1);
    for z=1:N
        tot_var(z)=sqrt(sum(eig(vect2mat_manual(cov_mat(z,:),D))));%vec2mat
    end
    minTV       = min(tot_var);
    maxTV       = max(tot_var);
    [uniqTvar, ~, IDg]=unique(tot_var);
    numUNQ=numel(uniqTvar);
    %For each unc. group calculate inter dist table
    for i=numUNQ:-1:2
        tempIDs  = IDg==i;  % Source IDs
        tempIDt  = IDg<i;   % Target IDs
        numS     = sum(tempIDs); % Number of source points
        AtempIDs = find(tempIDs);
        AtempIDt = find(tempIDt);
        
        if(mod(numUNQ-i,round(numUNQ/100))==0)
            clc;
            disp(['Condensing: ' num2str(100-(100*i/numUNQ),'%.0f') '% done']);   
        end
        %disp(['GrpNO: ' num2str(i) ' NumEL: ' num2str(numS) ' : ' datestr(now)]);
        %% Define upper and lower bound as R(i)+R(i-1)
        vecWGH  = con_cat(tempIDs,D+1); 
        cvrTBL  = cov_mat(tempIDs,:);
        
        matS       = org_cat(tempIDs,:); % Source coordinates
        matT       = org_cat(tempIDt,:); % Target coordinates
        
        % Check all events within -4 +4 sigma of the source events
        matSub = matS+repmat(tot_var(tempIDs),[1 D])*4;
        matSlb = matS-repmat(tot_var(tempIDs),[1 D])*4;

        for k=1:numS % Loop through all sources
            matchID     = all([bsxfun(@lt,matT,matSub(k,:)) bsxfun(@gt,matT,matSlb(k,:))],2);
            AmatchID    = AtempIDt(matchID);
            numT        = sum(matchID); % Number of possible targets
            if(numT)
                %Sample each source event with a number of points;
                mS      = matS(k,:)';
                cvrS    = vect2mat_manual(cvrTBL(k,:),D);%vec2mat
                testS   = randnorm(numMO, mS, [], cvrS);
                per_pdf = zeros(numMO,numT+1);
                per_pdf(:,numT+1) = gaussianValue(testS, mS, cvrS)*vecWGH(k);
                mT      = cell(numT,1);
                cvrT    = cell(numT,1);
                for j=1:numT %Calculate likelihood for all target events
                   mT{j}        = org_cat(AmatchID(j),:)';
                   cvrT{j}      = reshape(cov_mat(AmatchID(j),:),D,D);
                   per_pdf(:,j) = gaussianValue(testS, mT{j}, cvrT{j})*wght_vec(AmatchID(j));
                end
                [~,chrome]  = max(per_pdf,[],2);
                num_e       = hist(chrome,1:numT+1);
                vecOVL      = (num_e'/numMO)*vecWGH(k);
                
                con_cat(AmatchID,D+1)   = con_cat(AmatchID,D+1) +...
                                                              vecOVL(1:numT);
                con_cat(AtempIDs(k),D+1)                   = vecOVL(numT+1);
            end
        end
    end
    
    if(DO_PLOT)
        figure;
        %Original Catalog
        subplot(1,2,1);
        if(D==2)
            scatter(org_cat(:,1),org_cat(:,2),20,tot_var,'o','filled');
            hold on;
            for i=1:N
                
                %Plot 2D ellipses
                clr = jet(100);
                cid =   1 + round(99*((tot_var(i)-minTV)/(maxTV-minTV)));
                h   =   plot_gaussian_ellipsoid(org_cat(i,:), vect2mat_manual(cov_mat(i,:),D),1);%vec2mat
                set(h,'Color',clr(cid,:));
         
            end
            tmpC=colorbar('SouthOutside');
            set(get(tmpC,'xlabel'),'String','Isotropic Variance');
            title([ 'Original Catalog: ' num2str(N)] );
            daspect([1 1 1]);
            axis tight;
            xL=get(gca,'xlim');
            yL=get(gca,'ylim');
            xL=xL+[-1 1]*max(tot_var);
            yL=yL+[-1 1]*max(tot_var);
            set(gca,'xlim',xL);
            set(gca,'ylim',yL);
            caxis([minTV maxTV]);
            %Condensed Catalog
            subplot(1,2,2);
            IDnz = con_cat(:,D+1)>0; %>0
            scatter(org_cat(IDnz,1),org_cat(IDnz,2),20,con_cat(IDnz,D+1),'o','filled');
            daspect([1 1 1]);
            set(gca,'xlim',xL);
            set(gca,'ylim',yL);
            colormap(jet);
            tmpC=colorbar('SouthOutside');
            set(get(tmpC,'xlabel'),'String','Weights');
            title(['Condensed Catalog: ' num2str(sum(IDnz)) '(w>0)']);
            caxis([0 prctile(con_cat(:,D+1),85)]);
        elseif(D==3)
            scatter3(org_cat(:,1),org_cat(:,2),org_cat(:,3),20,tot_var,'o','filled');
            hold on;
            %Plot 3D ellipsoids
            clr = jet(100);
            for i=1:N
                cid =   1 + round(99*((tot_var(i)-minTV)/(maxTV-minTV)));
                h   =   plot_gaussian_ellipsoid(org_cat(i,:), vect2mat_manual(cov_mat(i,:),D),1);%vec2mat
                set(h,'edgecolor','none','facealpha',0.3,'FaceColor',clr(cid,:));
            end
            tmpC=colorbar('SouthOutside');
            set(get(tmpC,'xlabel'),'String','Isotropic Variance');
            title([ 'Original Catalog: ' num2str(N)] );
            daspect([1 1 1]);
            axis tight;
            xL=get(gca,'xlim');
            yL=get(gca,'ylim');
            zL=get(gca,'zlim');
            caxis([minTV maxTV]);
            set(gca,'zdir','reverse')
            view(0,-90);
            %Condensed Catalog
            subplot(1,2,2);
            IDnz = con_cat(:,D+1)>0; %>0
            scatter3(org_cat(IDnz,1),org_cat(IDnz,2),org_cat(IDnz,3),20,con_cat(IDnz,D+1),'o','filled');
            daspect([1 1 1]);
            set(gca,'xlim',xL);
            set(gca,'ylim',yL);
            set(gca,'zlim',zL);
            colormap(jet);
            tmpC=colorbar('SouthOutside');
            set(get(tmpC,'xlabel'),'String','Weights');
            title(['Condensed Catalog: ' num2str(sum(IDnz)) '(w>0)']);
            caxis([0 prctile(con_cat(:,D+1),85)]);
            set(gca,'zdir','reverse')
            view(0,-90);
        end
    end
end


