function J_clus=calc_J_per_clus(xst,yst,xv_clus,yv_clus,Lst,vec_planest)

%fprintf('From rectdist %g %g\n',[k m]);

% k = index of hypocenter wanted
% m = index of fault plane wanted

% Construct vertices vector
v=[xv_clus' yv_clus'];

xb = mean(xst);
yb = mean(yst);

J_clus = 0;

for k=1:length(xst)
    
    % Choose hypocenter vector
    u =[xst(k) yst(k);xst(k) yst(k)];
    u1=[xst(k) yst(k)];
        
    % Construct distance vector from each vertex
    d=u-v;
    d1=[d(1,1) d(1,2)];
    d2=[d(2,1) d(2,2)];

    d1n=norm(d1);
    d2n=norm(d2);

    % Find all angles and perpendicular distances to edge
    % 12 Face
    %gam12=acos(d12./(d1n.*d2n))
    if abs(d1n - d2n)<1e-4 % hypocenter exactly at the center of the line.
        alph12=0;
        beta12=0;
        gam12=pi;
    else

        [alph12,beta12,gam12]=trisol(d2n,d1n,Lst,'r');
        a12=d1n.*d2n.*sin(gam12)./Lst;
    end

    % Find the perpendicular distance to the infinite plane
    vec(1:2)=vec_planest;
    u2=u1-[xb yb];
    dperp=abs(dot(u2,vec));

    % Find the position field for the hypocenter relative to the plane and find
    % the minimum distance.
    p2=pi/2.0;
    dmin=dperp;
    if (alph12 > p2); dmin=d1n; end
    if (beta12 > p2); dmin=d2n; end

    %  accumulate the global variance
    J_clus = J_clus + dmin.*dmin;
    
end

% variance for the cluster
J_clus = J_clus./round(length(xst));

return;
