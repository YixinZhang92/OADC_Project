function OADC_2D_on_proj_hypos()

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
global index use_glo_var con_tol Kfaults

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
PLOT_FLAG1 = 1;

%******************** Big Loop over Kfaults *******************************
while Kfaults <= kmax
    %fprintf('*********************************************************\n\n');
    %fprintf('Kfaults= %i\n\n',Kfaults);

    %  form clusters of seismicity using present number of random faults.
    %  Much of the work is done here
    %fprintf('** Forming clusters of seismicity using present number of random faults **\n\n');

    if use_glo_var == 1
        JFINAL=faultcluster_2D(con_tol,Kfaults);
    else
        JFINAL=Copy_of_faultcluster_2D(con_tol,Kfaults);
    end

    %Kfaults;
    %Nt(1:Kfaults);

    %  plot initial planes
            if PLOT_FLAG1 == 1
                picname=strcat('Iteration',num2str(Kfaults),'Model');
                %datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
                datplot_2D(xs,ys,Kfaults,xv,yv,picname,simul_tag);
            end




    %  test to see if fit is within the error.  Look at largest lambda3
    %  eigenvalue
    lambda3max=max(lambda3);

    if lambda3max <= err_av

        %  print the good news
        %fprintf('Fault model converged to within error!\n');
        SOL_FLAG=1;
        Kfaults_good=Kfaults;
        Kfaults=kmax+1;

    else
        % split the thickest fault into two new random fault planes
        if Kfaults < kmax          
            %fprintf('*********************************************************\n\n');
            %fprintf('Splitting thickest fault, Kfaults= %i + 1\n',Kfaults);

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

            %textprogressbar('Determining the best random fault configurations: ');  

            splitfault_Nloop_2D

            %textprogressbar('done');

            JFINAL;

            % Determine the iteration with the minimum value of JFINAL, and
            % store its corresponding arrays
            [~,index] = sort(JFINAL);

            %  print the good news
%                     fprintf('Lambda3 of all configs: ');
%                     fprintf('[');fprintf('%8.4f', JFINAL);fprintf(']')
%                     fprintf('\n')
%                     fprintf('Best configuration found!\n');

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

            
            
%             
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
    %fprintf('Analysis Complete, fault model is not optimal\n');
    Kfaults=Kfaults-1;
else
    Kfaults=Kfaults_good;
end

% % saving all variables to file
% savevar_filename = [simul_tag '.saved_variables.mat'];
% save(savevar_filename)

%Strike
% 
% %  plot final planes
% picname='Final Model';
% %datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
% datplot_2D(xs,ys,Kfaults,xv,yv,picname,simul_tag); 

