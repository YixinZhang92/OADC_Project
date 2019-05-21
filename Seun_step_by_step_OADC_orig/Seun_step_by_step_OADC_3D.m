%function Seun_step_by_step_OADC_3D()
% Performing OADC_3D step by step by Seun using Chuck's functions
clear all; close all; clc


%***************** Read Catalog of Hypocenters ****************************
% data = textread('src.10');
% xs = data(:,5);
% ys = data(:,6);
% zs = -data(:,7);


%% OPEN HYPOCENTER FILE AND PLOT
%infile = 'result_declustered_hypo_0.8.txt';
%infile = 'result_collapsed_hypo_0.8.txt';
infile = 'result_declustered_collapsed_hypo_0.8.txt';

fid=fopen(infile,'r');
[data,count]=fscanf(fid,'%g %g %g',[3,inf]);
fclose(fid);

data=data'; [N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

% figure
% scatter3(xs,ys,zs,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); 
%(Coordiante (0,0) is at the lower right corner)
% grid MINOR; hold on




%% OADC
[xv,yv,zv,vec_plane,L,W] = Copy_of_randfaults(xs,ys,zs);

%  plot initial planes
picname='Initial Model';
datplot(xs,ys,zs,1,xv,yv,zv,picname);


%% ******************** Big Loop over Kfaults *******************************
%  form clusters of seismicity using present number of random faults.
%  Much of the work is done here
con_tol=0.1;  %  units usually in km
Kfaults = 1;

N = length(xs);

for k=1:N   %per hypocenter
        
    for m=1:Kfaults  %per fault plane
            
        dst(m)=Copy_of_rectdist(xs(k),ys(k),zs(k),xv,yv,zv,k,m,L,W);
            
    end
        
    dst;
    
    %  find the closest fault plane
    [val,index]=min(dst);
        
    Nt(index)=Nt(index) + 1;
    xt(index,Nt(index))=xs(k);
    yt(index,Nt(index))=ys(k);
    zt(index,Nt(index))=zs(k);
    
    %  accumulate the global variance
    J=J+val.*val;
        
end;








JFINAL=faultcluster(con_tol,Kfaults);





