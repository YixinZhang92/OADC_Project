%function OADC_3D_Seun_editted_7_2_2019()
%function OADC_3D_Seun_editted(kmin,kmax,err_av,infile)

%  Implementation of 3-D Optimal Anisotropic Dynamic Clustering from

%  Ouillon, Ducorbier, and Sornette (2008).  Automatic reconstruction of
%  fault networks from seismicity catalogs: Three-dimensional optimal
%  anisotropic dynamic clustering, JGR, 113, B01306,
%  doi:10.1029/2007JB00503.

%  specify:

%       kmin = starting number of fault planes
%       kmax = maximum number of fault planes analyzed
%       err_av = average hypocentral error in km
%       infile = file containing (x,y,z) positions of hypocenters

% <command execution_time = "8810">OADC_3D(1,7,0.01,'testdata.txt')</command>
% <command execution_time = "4459">OADC_3D(1,3,0.5,'testdata.txt')</command>
% <command execution_time = "5385">OADC_3D(1,4,0.5,'testdata.txt')</command>
% OADC_3D_Seun_editted(1,1,0.1,'result_declustered_collapsed_hypo_0.8.txt')


close all; clear all; clc;

%rng shuffle
rng('shuffle');

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global Strike Dip

%close all;
%kmin = 1;kmax=3;err_av=0.5;infile='testdata.txt';

kmin=1;kmax=3;err_av=0.01;%0.001;
N_loop = 10;
infile= 'COLCUM.20F_hypos.txt';
%infile= 'testdata.txt';
%infile= 'declustered_collapsed_hypo_0.6.txt';
simul_tag = 'Simul.1';

% remove previous calculations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

%infile= 'testdata.txt';


%********************** Set Parameters ************************************
%   Fault length scale for random faults.  Will be between 0 and fscale in
%   km
fscale=50.0;
%   Convergence tolerance value for the clustering algorithm in
%   'faultcluster'.  Represents the smallest change in global variance with
%   hypocenter clustering iteration.  The clustering process will stop once
%   the change in global variance with iteration drops to this value or
%   smaller.
con_tol=0.001;  %  units usually in km
PLOT_FLAG1=1;   % =0, no intermediate loop plots of data and planes

%***************** Read Catalog of Hypocenters ****************************
read_catalog(infile);

%********************** Initialize Space **********************************
init_space(kmax);

%******************* Initialize random faults *****************************
FAULT_FLAG=0;   % Initialization, use all hypocenters
Copy_of_randfaults(kmin,FAULT_FLAG);
%randfaults(kmin,FAULT_FLAG);

%  plot initial planes
picname='Initial Model';
Copy_of_datplot(xs,ys,zs,kmin,xv,yv,zv,picname,simul_tag);

SOL_FLAG=0;
Kfaults=kmin;

%% ******************** Big Loop over Kfaults *******************************
while Kfaults <= kmax
    fprintf('*********************************************************\n\n');
    fprintf('Kfaults= %i\n\n',Kfaults);
    
    %  form clusters of seismicity using present number of random faults.
    %  Much of the work is done here
    fprintf('** Forming clusters of seismicity using present number of random faults **\n\n');
    [JFINAL,gvar]=Copy_of_faultcluster_using_lamda3(con_tol,Kfaults);
    
    fprintf('Largest \x03bb_3 = %7.4f\n\n',JFINAL);
    
    %  plot initial planes
    if PLOT_FLAG1 == 1
        picname=strcat('Iteration',num2str(Kfaults),'Model');
        Copy_of_datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
    end
  
    %  test to see if fit is within the error.  Look at largest lambda3
    %  eigenvalue
    lambda3max=max(lambda3);
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
            
%             Copy_of_splitfault(Kfaults)
%             %Copy_of_splitfault_N_times(con_tol, Kfaults)
% 
%             % increase the fault number
%             Kfaults=Kfaults+1; 
        
            % Seun put the N loop here. I think it is better here.
            % load up arrays with good cluster parameters    
            xb_tmp=xb; yb_tmp=yb; zb_tmp=zb; % Barycenters
            xv_tmp=xv; yv_tmp=yv; zv_tmp=zv; % fault plane vertices
            xt_tmp=xt; yt_tmp=yt; zt_tmp=zt; % hypocenter location in a cluster
            vec_plane_tmp=vec_plane; % eigenvector that describes each plane
            Nt_tmp=Nt; % number of events in each trial cluster
            lambda3_tmp=lambda3; % minimum eigenvalue
            L_tmp=L; W_tmp=W; Strike_tmp=Strike; Dip_tmp=Dip; % fault plane parameters

     
            % Initialize loop temporary arrays
            [m,n] = size(xb); xb_tmp_i = zeros(m,n,N_loop); 
            yb_tmp_i = zeros(m,n,N_loop); zb_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(xv); xv_tmp_i = zeros(m,n,N_loop); 
                yv_tmp_i = zeros(m,n,N_loop); zv_tmp_i = zeros(m,n,N_loop);
            [~,n] = size(xt); xt_tmp_i = zeros(Kfaults+1,n,N_loop); yt_tmp_i = ...
                zeros(Kfaults+1,n,N_loop); zt_tmp_i = zeros(Kfaults+1,n,N_loop);
            [m,n] = size(vec_plane); vec_plane_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(Nt); Nt_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(lambda3); lambda3_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(L); L_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(W); W_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(Strike); Strike_tmp_i = zeros(m,n,N_loop);
            [m,n] = size(Dip); Dip_tmp_i = zeros(m,n,N_loop);
                    
            textprogressbar('Determining the best random fault configurations: ');

            for i = 1: N_loop
             
                % Put back the temporary arrays into the original arrays
                % load up arrays with good cluster parameters     
                xb = xb_tmp; yb = yb_tmp; zb = zb_tmp; % Barycenters
                xv=xv_tmp; yv=yv_tmp; zv=zv_tmp; % fault plane vertices
                xt=xt_tmp; yt=yt_tmp; zt=zt_tmp; % hypocenter location in a cluster
                vec_plane=vec_plane_tmp; % eigenvector that describes each plane
                Nt=Nt_tmp; % number of events in each trial cluster
                lambda3=lambda3_tmp; % minimum eigenvalue
                L=L_tmp; W=W_tmp; Strike=Strike_tmp; Dip=Dip_tmp; % fault plane parameters
            
                Copy_of_splitfault(Kfaults)
                %Copy_of_splitfault_N_times(con_tol, Kfaults)

                % increase the fault number
                Kfaults=Kfaults+1;
                
                [JFINAL(i), gvar(i)]=Copy_of_faultcluster_using_lamda3(con_tol,Kfaults);
                
                Kfaults=Kfaults-1;
      
                % Store each parameter to know which to advance to the next
                % iteration.           
                % load up arrays with good cluster parameters    
                xb_tmp_i(:,:,i)=xb; yb_tmp_i(:,:,i)=yb; zb_tmp_i(:,:,i)=zb; % Barycenters
                xv_tmp_i(:,:,i)=xv; yv_tmp_i(:,:,i)=yv; zv_tmp_i(:,:,i)=zv; % fault plane vertices
                xt_tmp_i(:,:,i)=xt; yt_tmp_i(:,:,i)=yt; zt_tmp_i(:,:,i)=zt; % hypocenter location in a cluster
                vec_plane_tmp_i(:,:,i)=vec_plane; % eigenvector that describes each plane
                Nt_tmp_i(:,:,i)=Nt; % number of events in each trial cluster
                lambda3_tmp_i(:,:,i)=lambda3; % minimum eigenvalue
                L_tmp_i(:,:,i)=L; W_tmp_i(:,:,i)=W; Strike_tmp_i(:,:,i)=Strike;
                Dip_tmp_i(:,:,i)=Dip; % fault plane parameters

                perc = (i/N_loop)*100;
                textprogressbar(perc);
                pause(0.1);
            end
            
            textprogressbar('done');
        
            % Determine the iteration with the minimum value of JFINAL, and
            % store its corresponding arrays
            [~,index] = sort(JFINAL);
            
            JFINAL
            gvar
            
            %  print the good news
            fprintf('Lambda3 of all configs: ');
            fprintf('[');fprintf('%8.4f', JFINAL);fprintf(']')
            fprintf('\n')
            fprintf('Best configuration found!\n');

            % load up arrays with good cluster parameters    
            xb = xb_tmp_i(:,:,index(1)); yb = yb_tmp_i(:,:,index(1)); 
            zb=zb_tmp_i(:,:,index(1)); % Barycenters
            xv = xv_tmp_i(:,:,index(1)); yv=yv_tmp_i(:,:,index(1)); 
            zv=zv_tmp_i(:,:,index(1)); % fault plane vertices
            xt=xt_tmp_i(:,:,index(1)); yt=yt_tmp_i(:,:,index(1)); 
            zt=zt_tmp_i(:,:,index(1)); % hypocenter location in a cluster
            vec_plane=vec_plane_tmp_i(:,:,index(1)); % eigenvector that describes each plane
            Nt=Nt_tmp_i(:,:,index(1)); % number of events in each trial cluster
            lambda3=lambda3_tmp_i(:,:,index(1)); % minimum eigenvalue
            L=L_tmp_i(:,:,index(1)); W=W_tmp_i(:,:,index(1)); 
            Strike=Strike_tmp_i(:,:,index(1)); Dip=Dip_tmp_i(:,:,index(1)); % fault plane parameters
            
            % increase the fault number
            Kfaults=Kfaults+1;     
            
            % plot new planes with data
            if PLOT_FLAG1 == 1
                picname='Model from fault split';
                Copy_of_datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
            end
            
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

%outfaults;
%% saving all variables to file
savevar_filename = [simul_tag '.saved_variables.mat'];
save(savevar_filename)

%%  plot final planes
picname='Final Model';
Copy_of_datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);

picname='Final Model 2'; n = 2; xv_disp = xv_tmp_i(:,:,index(n));
yv_disp = yv_tmp_i(:,:,index(n)); zv_disp = zv_tmp_i(:,:,index(n));
Copy_of_datplot(xs,ys,zs,Kfaults,xv_disp,yv_disp,zv_disp,picname,simul_tag);

picname='Final Model 3'; n = 3; xv_disp = xv_tmp_i(:,:,index(n));
yv_disp = yv_tmp_i(:,:,index(n)); zv_disp = zv_tmp_i(:,:,index(n));
Copy_of_datplot(xs,ys,zs,Kfaults,xv_disp,yv_disp,zv_disp,picname,simul_tag);

picname='Final Model 4'; n = 4; xv_disp = xv_tmp_i(:,:,index(n));
yv_disp = yv_tmp_i(:,:,index(n)); zv_disp = zv_tmp_i(:,:,index(n));
Copy_of_datplot(xs,ys,zs,Kfaults,xv_disp,yv_disp,zv_disp,picname,simul_tag);

picname='Final Model 5'; n = 5; xv_disp = xv_tmp_i(:,:,index(n));
yv_disp = yv_tmp_i(:,:,index(n)); zv_disp = zv_tmp_i(:,:,index(n));
Copy_of_datplot(xs,ys,zs,Kfaults,xv_disp,yv_disp,zv_disp,picname,simul_tag);

picname='Final Model 6'; n = 6; xv_disp = xv_tmp_i(:,:,index(n));
yv_disp = yv_tmp_i(:,:,index(n)); zv_disp = zv_tmp_i(:,:,index(n));
Copy_of_datplot(xs,ys,zs,Kfaults,xv_disp,yv_disp,zv_disp,picname,simul_tag);

%% Move all figures to a folder name with the simul_tag
eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',simul_tag, '*',simul_tag,'_results'))

