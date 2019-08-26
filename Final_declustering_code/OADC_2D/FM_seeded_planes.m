function FM_seeded_planes(n0,FAULT_FLAG)
%  Construct n0 faults using nearby focal mechanisms when we are plitting
%  the thickest fault in OADC_3D
%  just need the vertice locations

%  FAULT_FLAG = kthick (from splitfault.m),  Split the thickest fault into
%               two random planes of length L/2

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster


global xc yc vec_plane xb_old yb_old xs ys N Nc
global L xv yv L_old xv_old yv_old fscale
global xt yt Nt xb yb lambda3
global Strike FM_file dist2FM_threshold

% Load FM data set
fid=fopen(FM_file,'r');
% Check if file exists. If not run the random fault generator function
if fid == -1
    randfaults(n0,FAULT_FLAG) 
else
    
    [data,~]=fscanf(fid,'%g %g %g',[6,inf]); 
    fclose(fid);

    data=data'; [N_fm,~]=size(data);

    %  FM hypocentral locations, N=number of hypocenters
    xsf=data(1:N_fm,1);
    ysf=data(1:N_fm,2);
    zsf=data(1:N_fm,3);
    strikes=data(1:N_fm,4);
    dips = data(1:N_fm,5);
    rakes = data(1:N_fm,6);

    xthick = xt(FAULT_FLAG,1:Nt(FAULT_FLAG));
    ythick = yt(FAULT_FLAG,1:Nt(FAULT_FLAG));
    zthick = zt(FAULT_FLAG,1:Nt(FAULT_FLAG));

    % Find the closet FM to the earthquake in the thick cluster
    kk= 0; dst = zeros(1,N_fm);
    for i = 1:Nt(FAULT_FLAG) % per hypocenter

        for m=1:N_fm  %per FM
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
    % fprintf('No of FM near thick cluster = %i\n\n',kk);
    if kk==0 % No FM near the thick fault
        randfaults(n0,FAULT_FLAG)   
    else
        randfaults_using_FM(n0,FAULT_FLAG,strikes_fm,dips_fm,kk,xs_fm,ys_fm,zs_fm);
    end

end
