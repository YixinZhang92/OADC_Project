close all; clc; clear all;
tic

global xc yc zc vec_plane xb_old yb_old zb_old orig_xs orig_ys orig_zs N Nc xs ys
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global Strike Dip FM_file dist2FM_threshold dip_threshold
global xb_tmp_i yb_tmp_i zb_tmp_i
global xv_tmp_i yv_tmp_i zv_tmp_i
global xt_tmp_i yt_tmp_i zt_tmp_i
global vec_plane_tmp_i 
global Nt_tmp_i lambda3_tmp_i
global L_tmp_i W_tmp_i Strike_tmp_i Dip_tmp_i
global index use_glo_var con_tol Kfaults

kmin = 1; kmax=7; err_av=0.2;
N_loop = 1; simul_tag = 'Simul.manyFaults.OADC'; use_glo_var = 2;
FM_file='FM_dataset.csv'; dist2FM_threshold = 1; dip_threshold = 0;
%infile = 'Simul.1_hypos.txt';
infile = 'testdata.txt';
%infile = 'cluster3.txt';
%infile = 'Simul.now_ALL_hypos_hypos.txt';
%infile = 'Simul.2_ALL_hypos_hypos.txt';
%infile = 'CSZ_hypos.txt'; 
%infile='COLCUM.20F_hypos.txt';

rng('shuffle');

% remove previous simulations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

% Save diary
diaryname = [simul_tag '.myDiaryFile.txt'];
diary(diaryname) 

%% ********************** Set Parameters ************************************
%   Fault length scale for random faults. Will be between 0 and fscale in km
fscale=50.0;
%   Convergence tolerance value for the clustering algorithm in
%   'faultcluster'.  Represents the smallest change in global variance with
%   hypocenter clustering iteration.  The clustering process will stop once
%   the change in global variance with iteration drops to this value or
%   smaller.
con_tol=0.01;  %  units usually in km
PLOT_FLAG1=1; % =0, no intermediate loop plots of data and planes

ncount=0;
for az=0:20:179
    for el=-90:20:90
        ncount = ncount+1;
        %***************** Read Catalog of Hypocenters ****************************
        read_catalog_P3D(infile,simul_tag,0);

        % project the 3D hypos on 2D planes specified by azimuth and
        % elevations of viewpoint
        [xs,ys] = proj_of_3D_to_2D_plane(orig_xs,orig_ys,orig_zs,az,el);
        
%         X = [xs' ys'];
%         plot(xs,ys,'o');grid MINOR; shg
%         xlim([-5 5]);
%         ylim([-10 10]);
%         pause(1)       
        
        %********************** Initialize Space **********************************
        init_space_2D(kmax);

        %******************* Initialize random faults *****************************
        FAULT_FLAG=0;   % Initialization, use all hypocenters

        %initial_faults_using_kmeans_2D(kmin);
        randfaults_2D(kmin,FAULT_FLAG);

%         %  plot initial planes
%         picname='Initial Model';
%         datplot_2D(xs,ys,kmin,xv,yv,picname,simul_tag);

        SOL_FLAG=0;
        Kfaults=kmin;
        

        %******************** Big Loop over Kfaults *******************************
        while Kfaults <= kmax
            fprintf('*********************************************************\n\n');
            fprintf('Kfaults= %i\n\n',Kfaults);

            %  form clusters of seismicity using present number of random faults.
            %  Much of the work is done here
            fprintf('** Forming clusters of seismicity using present number of random faults **\n\n');

            if use_glo_var == 1
                JFINAL=faultcluster_2D(con_tol,Kfaults)
            else
                JFINAL=Copy_of_faultcluster_2D(con_tol,Kfaults)
            end
            
            Kfaults;
            Nt(1:Kfaults);

            %  plot initial planes
%             if PLOT_FLAG1 == 1
%                 picname=strcat('Iteration',num2str(Kfaults),'Model');
%                 %datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
%                 datplot_2D(xs,ys,Kfaults,xv,yv,picname,simul_tag);
%             end

            %  test to see if fit is within the error.  Look at largest lambda3
            %  eigenvalue
            lambda3max=max(lambda3)
            if lambda3max <= err_av

                %  print the good news
                fprintf('Fault model converged to within error!\n');
                SOL_FLAG=1;
                Kfaults_good=Kfaults;
                Kfaults=kmax+1;

            else
                % split the thickest fault into two new random fault planes
                if Kfaults < kmax          
                    fprintf('*********************************************************\n\n');
                    fprintf('Splitting thickest fault, Kfaults= %i + 1\n',Kfaults);

                    %             splitfault(Kfaults);
                    %          ..

                    %             % increase the fault number
                    %             Kfaults=Kfaults+1;

                    % OADC_3D will split the thickest fault N_loop times to find the
                    % configuration with the best fit.           
                    % load up arrays with good cluster parameters    
                    xb_tmp=xb; yb_tmp=yb; %zb_tmp=zb; % Barycenters
                    xv_tmp=xv; yv_tmp=yv;% zv_tmp=zv; % fault plane vertices
                    xt_tmp=xt; yt_tmp=yt; %zt_tmp=zt; % hypocenter location in a cluster
                    vec_plane_tmp=vec_plane; % eigenvector that describes each plane
                    Nt_tmp=Nt; % number of events in each trial cluster
                    lambda3_tmp=lambda3; % minimum eigenvalue
                    L_tmp=L; %W_tmp=W; 
                    Strike_tmp=Strike; %Dip_tmp=Dip; % fault plane parameters

                    % Initialize loop temporary arrays
                    [m,n] = size(xb); xb_tmp_i = zeros(m,n,N_loop); 
                    yb_tmp_i = zeros(m,n,N_loop); %zb_tmp_i = zeros(m,n,N_loop);
                    [m,n] = size(xv); xv_tmp_i = zeros(m,n,N_loop); 
                    yv_tmp_i = zeros(m,n,N_loop); %zv_tmp_i = zeros(m,n,N_loop);
                    [~,n] = size(xt); xt_tmp_i = zeros(Kfaults+1,n,N_loop); yt_tmp_i = ...
                    zeros(Kfaults+1,n,N_loop); %zt_tmp_i = zeros(Kfaults+1,n,N_loop);
                    [m,n] = size(vec_plane); vec_plane_tmp_i = zeros(m,n,N_loop);
                    [m,n] = size(Nt); Nt_tmp_i = zeros(m,n,N_loop);
                    [m,n] = size(lambda3); lambda3_tmp_i = zeros(m,n,N_loop);
                    [m,n] = size(L); L_tmp_i = zeros(m,n,N_loop);
                    %[m,n] = size(W); W_tmp_i = zeros(m,n,N_loop);
                    [m,n] = size(Strike); Strike_tmp_i = zeros(m,n,N_loop);
                    %[m,n] = size(Dip); Dip_tmp_i = zeros(m,n,N_loop);

                    textprogressbar('Determining the best random fault configurations: ');  

                    splitfault_Nloop_2D

                    textprogressbar('done');

                    JFINAL;

                    % Determine the iteration with the minimum value of JFINAL, and
                    % store its corresponding arrays
                    [~,index] = sort(JFINAL);

                    %  print the good news
                    fprintf('Lambda3 of all configs: ');
                    fprintf('[');fprintf('%8.4f', JFINAL);fprintf(']')
                    fprintf('\n')
                    fprintf('Best configuration found!\n');

                    % load up arrays with good cluster parameters    
                    xb = xb_tmp_i(:,:,index(1)); yb = yb_tmp_i(:,:,index(1)); 
                    %zb=zb_tmp_i(:,:,index(1)); % Barycenters
                    xv = xv_tmp_i(:,:,index(1)); yv=yv_tmp_i(:,:,index(1)); 
                    %zv=zv_tmp_i(:,:,index(1)); % fault plane vertices
                    xt=xt_tmp_i(:,:,index(1)); yt=yt_tmp_i(:,:,index(1)); 
                    %zt=zt_tmp_i(:,:,index(1)); % hypocenter location in a cluster
                    vec_plane=vec_plane_tmp_i(:,:,index(1)); % eigenvector that describes each plane
                    Nt=Nt_tmp_i(:,:,index(1)); % number of events in each trial cluster
                    lambda3=lambda3_tmp_i(:,:,index(1)); % minimum eigenvalue
                    L=L_tmp_i(:,:,index(1));% W=W_tmp_i(:,:,index(1)); 
                    Strike=Strike_tmp_i(:,:,index(1));% Dip=Dip_tmp_i(:,:,index(1)); % fault plane parameters

                    % increase the fault number
                    Kfaults=Kfaults+1;     

%                     %  plot new planes with data
%                     if PLOT_FLAG1 == 1
%                         picname=['Model from fault split ' num2str(Kfaults)];
%                         %datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
%                         datplot_2D(xs,ys,Kfaults,xv,yv,picname,simul_tag);
%                     end

                else
                    Kfaults=Kfaults+1;

                end
            end
        end

        %  Output the final fault model (knowing it is not optimal)
        if SOL_FLAG == 0
            fprintf('Analysis Complete, fault model is not optimal\n');
            Kfaults=Kfaults-1;
        else
            Kfaults=Kfaults_good;
        end

        % saving all variables to file
        savevar_filename = [simul_tag '.saved_variables.mat'];
        save(savevar_filename)

        Strike

%         %  plot final planes
%         picname='Final Model';
%         %datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
%         datplot_2D(xs,ys,Kfaults,xv,yv,picname,simul_tag); 
   
    end
end

toc





















%         
%         idx = rev_dbscan(X,0.5,10); % The default distance metric is Euclidean distance
%         %idx = DBSCAN(X,1,10); % The default distance metric is Euclidean distance
% 
%         %gscatter(X(:,1),X(:,2),idx); 
%         %title('DBSCAN Using Euclidean Distance Metric')
% 
%         %figure;
%         colors = {'r.', 'b.', 'g.', 'k.', 'y.', 'b.', 'k.', 'g.', 'c.', 'm.', ...
%                   'r.', 'b.', 'g.', 'k.', 'b.', 'k.', 'g.', 'c.', 'm.', };
% 
%         n0 = max(idx);
%         kmin = max(idx);
% 
%         for i=1:n0
%             %plot(X(idx==i,1),X(idx==i,2),colors{mod(i,11)+1},'MarkerSize',12); hold on
% 
%             % Calculate new barycenters for each cluster
%             xb(i) = mean(X(idx==i,1));
%             yb(i) = mean(X(idx==i,2));
%             Nt(i)= length(X(idx==i,1));
% 
%         end
%         %plot(X(idx<0,1),X(idx<0,2),'r.','MarkerSize',12); hold on
%         %title('DBSCAN Using Euclidean Distance Metric')
