% Declustering Algorithms driver
clear all; close all; clc; tic;

global xs ys zs

%infile = 'Simul.now_ALL_hypos_hypos.txt';
%infile = 'Simul.1_hypos.txt';
infile = 'CSZ_hypos.txt';
%infile = 'testdata.txt';
    
min_hypo_per_cluster = 3; 
err_x= 0.15; err_y= 0.15; err_z= 0.35;
simul_tag = 'Cum.with.OADC3D';

plot_FM = 0; % 0 - No FM; 1- FM on fault model only; 2- FM on all plots
infile_FM = 'Combined_Dataset_MT_PL.csv'; quality_threshold = 3; 
use_mag_size=1; FM_size_or_factor = 0.5;


%% Begin Clustering

% calculate volume threshold 
%vol_thresh = 4*sqrt(((3*err)^6)/16); %
vol_thresh = 4*(4/3)*pi*err_x*err_y*err_z; %
%vol_thresh = 0.2;

% remove previous calculations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

% Figure 1
read_catalog(infile,simul_tag)

if plot_FM == 2
    hold on;
    plot_FM_on_fault_model(infile_FM,quality_threshold,use_mag_size,FM_size_or_factor);
end

N = length(xs); clus = 0;
v(1:N)=0.0; index(1:N)=0.0; 

for k=1:N
    
    x0=xs(k); y0=ys(k); z0=zs(k);
    
    % calculate distance from reference source location
    dist=sqrt((xs-x0).^2 + (ys-y0).^2 + (zs-z0).^2);
    
    % sort by ascending distance
    [~,I]=sort(dist);
    
    %  use first 4 locations to compute volume of the tetrahedron
    %  associated with this source point
    x0=xs(I(1)); x1=xs(I(2)); x2=xs(I(3)); x3=xs(I(4));
    y0=ys(I(1)); y1=ys(I(2)); y2=ys(I(3)); y3=ys(I(4));
    z0=zs(I(1)); z1=zs(I(2)); z2=zs(I(3)); z3=zs(I(4));
        
    matr = [x0 y0 z0 1; x1 y1 z1 1; x2 y2 z2 1; x3 y3 z3 1];
    v(k) = abs(det(matr))./6;
        
    if v(k)<vol_thresh
        max_index = max([index(I(1)) index(I(2)) index(I(3)) index(I(4))]);
        
        if max_index == 0
            clus = clus + 1;
            
            index(I(1)) = clus;
            index(I(2)) = clus;
            index(I(3)) = clus;
            index(I(4)) = clus;
        else
            
            index(I(1)) = max_index;
            index(I(2)) = max_index;
            index(I(3)) = max_index;
            index(I(4)) = max_index;
        end 
    end
end

% remove hypocenters that have vol < vol_thresh and 
% clusters with less than min_hypo_per_cluster
index(v>vol_thresh)=0;
ind_clus = unique(index);

for i=2:length(ind_clus)
    ind = ind_clus(i);
    l_x = length(xs(index==ind));
    
    if l_x < min_hypo_per_cluster
        index(index==ind)=0;
    end
end

kk = 0;
for k=1:length(ind_clus)-1
    ind = ind_clus(k+1);  
    
    xst=xs(index==ind); yst=ys(index==ind); zst=zs(index==ind);

    % compute the covariance matrix for this cluster
    Cxy=cov([xst' yst' zst'],0);
    
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
    XX = [xst' yst' zst'];
    
    [Ltn(kk),Wtn(kk),Striket(kk),Dipt(kk),xvt(kk,:),yvt(kk,:),zvt(kk,:)] = ...
        fltplane(XX,V,D,xbt(kk),ybt(kk),zbt(kk));
            
    % save the plane unit normal vector and eigenvalue
    vec_plane_t(kk,1:3)=V(1:3,1);
    lambda3_t(kk)=sqrt(D(1,1));
end

%% Ploting figures
Fig1 = figure('Name','Fault Models with hypocenters',...
    'Position', get(0, 'Screensize'));

len=0;
for i=2:length(ind_clus)
    ind = ind_clus(i);
    plot3(xs(index==ind),ys(index==ind),zs(index==ind),'o'); hold on;
end

plot3(xs(index==ind_clus(1)),ys(index==ind_clus(1)),...
    zs(index==ind_clus(1)),'k.'); hold on;

for m=1:kk
    fill3(xvt(m,1:4),yvt(m,1:4),zvt(m,1:4),...
        'w','FaceAlpha',0.2,'FaceColor',[0.5 0.5 0.5]);
end

if plot_FM > 0
    hold on;
    plot_FM_on_fault_model(infile_FM,quality_threshold,use_mag_size, FM_size_or_factor);
end

grid on; axis equal;
title('Fault Model with isolated hypocenters');
set(gca, 'fontsize', 18); shg

%% Printing figure to file
fig_filename = [simul_tag '.faultmodel.png'];
F1    = getframe(Fig1);
imwrite(F1.cdata, fig_filename, 'png')
savefig(Fig1,[fig_filename(1:end-4) '.fig'])

% Move all figures to a folder name with the simul_tag
eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',simul_tag, '*',simul_tag,'_results'))

toc























% syn = 0; % syn = 0 - Real hypocenters
%          %     = 1 - Synthetic hypocenters
% if syn == 1
%     % synthetic hypocenters
%     mu= 0; sigma=1; nhypos = 1000;
%     outfile='rnd_hypos_3D.txt';
%     infile = outfile;
%     [data,R] = create_syn_rnd_hypo_3D(nhypos,mu,sigma,outfile);
% 
% else
%     % Real data
% 
% 
% end