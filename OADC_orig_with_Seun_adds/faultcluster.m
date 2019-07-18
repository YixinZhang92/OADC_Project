function [JFINAL]=faultcluster(con_tol,n0)
%  faultcluster - cluster the seismicity around the present set of faults.
%  Clusters will be adjusted until the global variance does not change
%  within the convergence tolerance parameter 'con_tol'
%  n0 is the number of faults presently being considered

%  JFINAL = global variance of the fit

% Seun added the following:
% 1. Maximum number of iteration in cae the iteration is endless.
% 2. Removed display of J at each iteration.

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale

DJ=2.*con_tol;
kj=0;
max_iter = 20;

%fprintf('from faultcluster\n');

% clustering loop
while (DJ > con_tol) && (kj < max_iter)
    kj=kj+1;
    %  form clusters around present number of fault planes and compute
    %  global variance J
    J=pcluster(n0)
    
    %fprintf('J= %g\n',J);
    
    if kj == 1
        JOLD=J;
    else
        JNEW=J;
        DJ=abs(JOLD-JNEW);
        JOLD=JNEW;
    end
    
    %  Compute Cxy matrix, perform principal components analysis, create
    %  new fault planes for each cluster
    recalcfault(n0);
    
end

JFINAL=JNEW;

return


