function OADC_3D(kmin, kmax,err_av,infile)

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
% OADC_3D(7,0.5,'testdata.txt')
% OADC_3D(1,7,0.1,'syn_hypo_all.txt')

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale

close all; 
%kmin = 3;kmax=3;err_av=0.5;infile='testdata.txt';



%********************** Set Parameters ************************************
%   Fault length scale for random faults.  Will be between 0 and fscale in
%   km
fscale=50.0;
%   Convergence tolerance value for the clustering algorithm in
%   'faultcluster'.  Represents the smallest change in global variance with
%   hypocenter clustering iteration.  The clustering process will stop once
%   the change in global variance with iteration drops to this value or
%   smaller.
con_tol=0.1;  %  units usually in km
PLOT_FLAG1=1;   % =0, no intermediate loop plots of data and planes

%***************** Read Catalog of Hypocenters ****************************
read_catalog(infile);

%********************** Initialize Space **********************************
init_space(kmax);

%******************* Initialize random faults *****************************
FAULT_FLAG=0;   % Initialization, use all hypocenters

% Using random faults as initial fault planes
randfaults(kmin,FAULT_FLAG);

% Using initial fault planes from the clustering analysis
% [kmin,analy, value_counts] = clustering_analysis(xs, ys, zs);
% 
% faults_from_clustering_analy(kmin,analy, value_counts);


%  plot initial planes
picname='Initial Model';
datplot(xs,ys,zs,kmin,xv,yv,zv,picname);

SOL_FLAG=0;
Kfaults=kmin;



%******************** Big Loop over Kfaults *******************************
while (Kfaults <= kmax) && (SOL_FLAG == 0)
    Kfaults
    %  form clusters of seismicity using present number of random faults.
    %  Much of the work is done here
    JFINAL=faultcluster(con_tol,Kfaults)
    
    %  plot initial planes
    if PLOT_FLAG1 == 1;
        picname=strcat('Iteration',num2str(Kfaults),'Model');
        datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
    end

    %  test to see if fit is within the error.  Look at largest lambda3
    %  eigenvalue
    lambda3max=max(lambda3)
    if lambda3max <= err_av;
        
        %  print the good news
        fprintf('Fault model converged to within error!\n');
        SOL_FLAG=1;
        Kfaults_good=Kfaults;
        Kfaults=kmax+1;
         
       
        
    else
        % split the thickest fault into two new random fault planes
        if Kfaults < kmax;
            fprintf('Splitting thickest fault, Kfaults= %i +1\n',Kfaults);
            splitfault(Kfaults);
         
            % increase the fault number
            Kfaults=Kfaults+1;
         
            %  plot new planes with data
            if PLOT_FLAG1 == 1;
                picname='Model from fault split';
                datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
            end
            
        else
            Kfaults=Kfaults+1;

        end;
        
    end;
     
end;

%  Output the final fault model (knowing it is not optimal)
if SOL_FLAG == 0; 
    fprintf('Analysis Complete, fault model is not optimal\n');
    Kfaults=Kfaults-1;
else
    Kfaults=Kfaults_good;
end

%outfaults;

%  plot final planes
picname='Final Model';
datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
        
end

