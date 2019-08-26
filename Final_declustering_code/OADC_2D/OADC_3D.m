%function OADC_3D(kmin,kmax,err_av,N_loop,infile,simul_tag,use_glo_var)
%,FM_file,dist2FM_threshold
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

% <command execution_time="8810">OADC_3D(1,7,0.01,'testdata.txt')</command>
% <command execution_time="4459">OADC_3D(1,3,0.5,'testdata.txt')</command>
% <command execution_time="5385">OADC_3D(1,4,0.5,'testdata.txt')</command>

% Modified by Oluwaseun Idowu Fadugba (7/11/2019)
% 1.  Shuffled random seed generator at the start of the code.
% 2.  Number of times OADC_3D splits the thichest faults
% 3.  Simulation tag to copy all generated files and figures to a folder. Delete
%     folder with the same simulation tag i.e., previous runs with the same tag.
% 4.  Edited randfaults.m
% 5.  Edited read_catalog.m, datplot.m to included simul_tag, and save figures with 
%     the tag as the filename. He also included figure title and adjusted fontsize.
% 6.  Added maximum number of iteration in faultcluster.m in case the iteration is endless.
% 7.  Initialized xt, yt and zt variables in pcluster.m because it later becomes 
%     inconsistent with the total number of earthquakes in the catalog.       
% 8.  Changed the lambda3(k)=sqrt(12.*D(1,1)) to lambda3(k)=sqrt(D(1,1)) in recalcfault.m
% 9.  OADC_3D will split the thickest fault N_loop times to find the
%     configuration with the best fit.
% 10. OADC_3D now saves all variables to file in .mat format.
% 11. OADC_3D now prints the best 6 fault geometries, and print them to file.
% 12. Added FM_file, dist2FM_threshold, strike and dip, and some temporary variables 
%     associated with fault split to global variables in OADC_3D.m
% 13. Modified display in OADC_3D to show status bar.
% 14. OADC_3D can now find a stable fault geometry using lambda_3 in
%     addition to global variance, specify by the use_glo_var.
%     use_glo_var = 1; use global variance
%     use_glo_var = 2; use lambda_3
% 15. OADC_3D can now split the thickest fault using the focal mechanisms
%     of earthquakes, if given, instead of randomly-seeded planes.
% 16. Changed the equations for determining strike and dip accordingly for
%     strike to go from the north clockwisely.
% 17. Fixed a bug in dividing thick fault into two in randfaults_using_FM.m by moving the 
%     L/2 out of the k=1:n0 loop, and changed L2 to L22 cos there is another L2 later.
% 18. Checked if Cxy contains NaN. This happens if no hypocenter is close to
%     one of the thickest faults.
% 19  Save diary to file using "diary [simul_tag '.myDiaryFile.txt']"
% 20  Plot best 6 fault models and display their fault parameters.
% 21  OADC_3D do not allow faults without an earthquake.
% 22  Splitfault is now within a while loop, until we get N_loop geometries
%     with faults having dips greater than the specified dip_threshold.
% 23  Seun added hypos in cluster to the input parameters to fltplane.m
% 24  Seun changed the way to determine L and W of fault plane


% Try to resolve it in randfaults.m as well.
%
%
% Run by typing:
%               OADC_3D(1,3,0.5,10,'testdata.txt','Simul.1',2)
%               OADC_3D(1,2,1,10,'COLCUM.20F_hypos.txt','Test.1',2)
%
%               OADC_3D(1,2,1,10,'COLCUM.20F_hypos.txt','Test.1',2,'FM_dataset.csv',1)

close all; clc; clear all;
tic

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
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

kmin = 1; kmax=3; err_av=1;
N_loop = 6; simul_tag = 'Simul.5Faults.OADC'; use_glo_var = 1;
FM_file='FM_dataset.csv'; dist2FM_threshold = 1; dip_threshold = 10;
infile = 'Simul.1_hypos.txt';
%infile = 'Simul.2_ALL_hypos_hypos.txt';
%infile = 'CSZ_hypos.txt'; c
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

%***************** Read Catalog of Hypocenters ****************************
read_catalog(infile,simul_tag);

%********************** Initialize Space **********************************
init_space(kmax);

%******************* Initialize random faults *****************************
FAULT_FLAG=0;   % Initialization, use all hypocenters
randfaults(kmin,FAULT_FLAG);

%  plot initial planes
picname='Initial Model';
datplot(xs,ys,zs,kmin,xv,yv,zv,picname,simul_tag);

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
        
        JFINAL=faultcluster(con_tol,Kfaults)
    else
        
        JFINAL=Copy_of_faultcluster(con_tol,Kfaults)
    end
    Kfaults;
     Nt(1:Kfaults);

    %  plot initial planes
    if PLOT_FLAG1 == 1
        picname=strcat('Iteration',num2str(Kfaults),'Model');
        datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
    end

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

            
            splitfault_Nloop
            



%              for i = 1: N_loop       
%                 % Put back the temporary arrays into the original arrays
%                 % load up arrays with good cluster parameters     
%                 xb = xb_tmp; yb = yb_tmp; zb = zb_tmp; % Barycenters
%                 xv=xv_tmp; yv=yv_tmp; zv=zv_tmp; % fault plane vertices
%                 xt=xt_tmp; yt=yt_tmp; zt=zt_tmp; % hypocenter location in a cluster
%                 vec_plane=vec_plane_tmp; % eigenvector that describes each plane
%                 Nt=Nt_tmp; % number of events in each trial cluster
%                 lambda3=lambda3_tmp; % minimum eigenvalue
%                 L=L_tmp; W=W_tmp; Strike=Strike_tmp; Dip=Dip_tmp; % fault plane parameters
%                             
%                 splitfault(Kfaults)
%                          
%                 % Store each parameter to know which to advance to the next
%                 % iteration.           
%                 % load up arrays with good cluster parameters    
%                 xb_tmp_i(:,:,i)=xb; yb_tmp_i(:,:,i)=yb; zb_tmp_i(:,:,i)=zb; % Barycenters
%                 xv_tmp_i(:,:,i)=xv; yv_tmp_i(:,:,i)=yv; zv_tmp_i(:,:,i)=zv; % fault plane vertices
%                % xt_tmp_i(:,:,i)=xt; yt_tmp_i(:,:,i)=yt; zt_tmp_i(:,:,i)=zt; 
%                % hypocenter location in a cluster
%                 vec_plane_tmp_i(:,:,i)=vec_plane; % eigenvector that describes each plane
%                 Nt_tmp_i(:,:,i)=Nt; % number of events in each trial cluster
%                 lambda3_tmp_i(:,:,i)=lambda3; % minimum eigenvalue
%                 L_tmp_i(:,:,i)=L; W_tmp_i(:,:,i)=W; Strike_tmp_i(:,:,i)=Strike;
%                 Dip_tmp_i(:,:,i)=Dip; % fault plane parameters
%   
%                 % increase the fault number
%                 Kfaults=Kfaults+1;
%     
%                 if use_glo_var == 1
%                     
%                     JFINAL(i)=faultcluster(con_tol,Kfaults);
%                 else
%                     
%                     JFINAL(i)=Copy_of_faultcluster(con_tol,Kfaults);
%                 end
%    
%                 % Reduce the number of fault because faultcluster.m have increased it by 1.
%                 Kfaults=Kfaults-1; 
%                  
% %                 % Store each parameter to know which to advance to the next
% %                 % iteration.           
% %                 % load up arrays with good cluster parameters    
% %                 xb_tmp_i(:,:,i)=xb; yb_tmp_i(:,:,i)=yb; zb_tmp_i(:,:,i)=zb; % Barycenters
% %                 xv_tmp_i(:,:,i)=xv; yv_tmp_i(:,:,i)=yv; zv_tmp_i(:,:,i)=zv; % fault plane vertices
% %                 xt_tmp_i(:,:,i)=xt; yt_tmp_i(:,:,i)=yt; zt_tmp_i(:,:,i)=zt; 
% %                   hypocenter location in a cluster
% %                 vec_plane_tmp_i(:,:,i)=vec_plane; % eigenvector that describes each plane
% %                 Nt_tmp_i(:,:,i)=Nt; % number of events in each trial cluster
% %                 lambda3_tmp_i(:,:,i)=lambda3; % minimum eigenvalue
% %                 L_tmp_i(:,:,i)=L; W_tmp_i(:,:,i)=W; Strike_tmp_i(:,:,i)=Strike;
%                  %Dip_tmp_i(:,:,i)=Dip; % fault plane parameters
%                  
%                  Dip
% min(Dip_tmp_i(:,:,i))
% 
%                 perc = (i/N_loop)*100;
%                 textprogressbar(perc);
%                 pause(0.1);
%                 
%             end
           







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
            
            %  plot new planes with data
            if PLOT_FLAG1 == 1
                picname=['Model from fault split ' num2str(Kfaults)];
                datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
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

% saving all variables to file
savevar_filename = [simul_tag '.saved_variables.mat'];
save(savevar_filename)

Strike
Dip

%%  plot final planes
picname='Final Model';
datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);

fprintf('*********************************************************\n\n');
fprintf('Final Model results 1:\n\n');
plot_figures_for_OADC(1)

fprintf('*********************************************************\n\n');
fprintf('Final Model results 2:\n\n');
plot_figures_for_OADC(2)

fprintf('*********************************************************\n\n');
fprintf('Final Model results 3:\n\n');
plot_figures_for_OADC(3)

fprintf('*********************************************************\n\n');
fprintf('Final Model results 4:\n\n');
plot_figures_for_OADC(4)

fprintf('*********************************************************\n\n');
fprintf('Final Model results 5:\n\n');
plot_figures_for_OADC(5)

fprintf('*********************************************************\n\n');
fprintf('Final Model results 6:\n\n');
plot_figures_for_OADC(6)

%% Move all figures to a folder name with the simul_tag
eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',simul_tag, '*',simul_tag,'_results'))

toc
%end