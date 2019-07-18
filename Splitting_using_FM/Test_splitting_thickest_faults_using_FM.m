function Test_splitting_thickest_faults_using_FM()
clear all; close all; clc;
load('Test.1.saved_variables.mat');

FM_file = 'FM_dataset.csv';
FAULT_FLAG = 2;
n0 = 2;
addpath('/Users/oluwaseunfadugba/Documents/OADC_project/OADC_orig_with_Seun_adds');
dist2FM_threshold = 1;
kmin=2;

%% Load FM data set
fid=fopen(FM_file,'r');
% Check if file exists. If not run the random fault generator function
if fid == -1
    randfaults(n0,FAULT_FLAG) 
else
    
    [data,~]=fscanf(fid,'%g %g %g',[6,inf]); 
    fclose(fid);

    data=data'; [N,~]=size(data);

    %  FM hypocentral locations, N=number of hypocenters
    xsf=data(1:N,1);
    ysf=data(1:N,2);
    zsf=data(1:N,3);
    strikes=data(1:N,4);
    dips = data(1:N,5);
    rakes = data(1:N,6);

    xthick = xt(FAULT_FLAG,1:Nt(FAULT_FLAG));
    ythick = yt(FAULT_FLAG,1:Nt(FAULT_FLAG));
    zthick = zt(FAULT_FLAG,1:Nt(FAULT_FLAG));

    % Find the closet FM to the earthquake in the thick cluster
    kk= 0; dst = zeros(1,N);
    for i = 1:Nt(FAULT_FLAG) % per hypocenter

        for m=1:N  %per FM
            dst(m)= sqrt((xsf(m)-xthick(i))^2 + (ysf(m)-ythick(i))^2 + (zsf(m)-zthick(i))^2);           
        end

        %  find the closest fault plane
        [mindist, index] = min(dst);

        if mindist < dist2FM_threshold 
            kk = kk + 1;

            xs_fm(kk) = xthick(i);
            ys_fm(kk) = ythick(i);
            zs_fm(kk) = zthick(i);

            strikes_fm(kk) = strikes(index);
            dips_fm(kk) = dips(index);
            rakes_fm(kk) = rakes(index);        
        end
    end

    % At this point, we know if there are FMs close to the thick cluster.
    % If no FM is present, OADC_3D will use random-seeded faults.
    if kk==0 % No FM near the thick fault
        randfaults(n0,FAULT_FLAG)   
    else
        randfaults_using_FM(n0,FAULT_FLAG,strikes_fm,dips_fm,kk,xs_fm,ys_fm,zs_fm);
    end

end


%  plot planes
picname='Initial Model';
datplot(xs,ys,zs,kmin,xv,yv,zv,picname,simul_tag);
shg
