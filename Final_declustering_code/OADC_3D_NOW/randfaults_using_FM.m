function randfaults_using_FM(n0,FAULT_FLAG,strikes_fm,dips_fm,kk,xs_fm,ys_fm,zs_fm)

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

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt yt zt Nt xb yb zb lambda3

%fprintf('From randfaults\n');

nbins= 20;
[freq,strike_freq] = hist(strikes_fm,nbins);
% figure
% hist(strikes_fm,nbins);

interval = (strike_freq(2)-strike_freq(1))/2;

[freqn,I]=sort(freq,'descend');
strike_freqn = strike_freq(I);

% Seun fixed a bug in dividing thick fault into two by moving the L/2 out
% of the k=1:n0 loop, and changed L2 to L22 cos there is another L2 later.
L22=L(FAULT_FLAG)./2.0;
W22=W(FAULT_FLAG);
    
for k = 1:n0

%     L(k)=W2 ;
%     W(k)=L2;
    
    %  make L2 >= W2

    if L22 <= W22
        L(k)=W22; 
        W(k)=L22;
    else
        L(k)=L22;
        W(k)=W22;
    end

    nb=randperm(kk);

    if kk == 1
        xb(k)=mean(xs_fm);
        yb(k)=mean(ys_fm);
        zb(k)=mean(zs_fm);
    else
        xb(k)=xs_fm(nb(1));
        yb(k)=ys_fm(nb(1));
        zb(k)=zs_fm(nb(1));
    end 

    % Use the dominant strike for the first fault and randomly chose the
    % other orientation.
%     if k == 1
%         strikes_fm_now = strike_freqn(k);
%         condition = strikes_fm>=(strikes_fm_now-interval) & ...
%             strikes_fm<=(strikes_fm_now+interval);
%         dips_fm_now = mean(dips_fm(condition));
% 
%     else
%         nb=randperm(length(freqn(freqn~=0)));
%         strikes_fm_now = strike_freqn(nb(1));
%         condition = strikes_fm>=(strikes_fm_now-interval) & ...
%             strikes_fm<=(strikes_fm_now+interval);
%         dips_fm_now = mean(dips_fm(condition));
%                                 
%     end

    % Use the dominant strike for the first fault and second dominant
    % strike for the second fault.
    strikes_fm_now = strike_freqn(k);
    condition = strikes_fm>=(strikes_fm_now-interval) & ...
        strikes_fm<=(strikes_fm_now+interval);
    dips_fm_now = mean(dips_fm(condition));
                        

    % Get unit vector from strike and dip.
    %  find the plane unit vector
    V(1) = -cosd(strikes_fm_now)*sind(dips_fm_now);
    V(2) =  sind(strikes_fm_now)*sind(dips_fm_now);
    V(3) = -cosd(dips_fm_now);
    vec_plane(k,1:3)= V;

    L2=L(k)./2.0;
    W2=W(k)./2.0;

    % compute vertice locations
    % First, for horizontal fault
    xvt = [-W2 -W2 W2 W2];
    yvt = [L2 -L2 -L2 L2];
    zvt = [0 0 0 0];

    % Fault corner coordinates have been created - Rotate the data into the
    % geographical coordinate system in order of rake, dip, strike
    % rotate into rake direction
    R=[xvt ; yvt ;zvt];
    [~,n] = size(R); rake =0;

    Drake=[cosd(rake) -sind(rake) 0 ;sind(rake) cosd(rake) 0 ;0 0 1];
    Rrake=Drake*R;

    xvt(1:n)=Rrake(1,1:n);
    yvt(1:n)=Rrake(2,1:n);
    zvt(1:n)=Rrake(3,1:n);

    % rotate into dip direction
    Ddip=[cosd(dips_fm_now) 0 sind(dips_fm_now); 0 1  0 ;...
         -sind(dips_fm_now) 0 cosd(dips_fm_now)];
    Rdip=Ddip*Rrake;

    xvt(1:n)=Rdip(1,1:n);
    yvt(1:n)=Rdip(2,1:n);
    zvt(1:n)=Rdip(3,1:n);

    % rotate into strike direction
    Dstrike=[cosd(strikes_fm_now) sind(strikes_fm_now) 0 ; ...
            -sind(strikes_fm_now) cosd(strikes_fm_now) 0 ; 0 0 1];
    Rstrike=Dstrike*Rdip;

    xvt(1:n)=Rstrike(1,1:n);
    yvt(1:n)=Rstrike(2,1:n);
    zvt(1:n)=Rstrike(3,1:n);

    xv(k,1:n)=xvt(1:n) + xb(k);
    yv(k,1:n)=yvt(1:n) + yb(k);
    zv(k,1:n)=zvt(1:n) + zb(k);
end