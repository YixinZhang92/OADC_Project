function [max_lamda3FINAL,J]=Copy_of_faultcluster_using_lamda3(con_tol,n0)
%  faultcluster - cluster the seismicity around the present set of faults.
%  Clusters will be adjusted until the global variance does not change
%  within the convergence tolerance parameter 'con_tol'
%  n0 is the number of faults presently being considered

%  JFINAL = global variance of the fit

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale

max_lamda3J=2.*con_tol;
kj=0;
max_iter = 20;

%fprintf('from faultcluster\n');

% clustering loop
while (max_lamda3J > con_tol) && (kj < max_iter)
    kj=kj+1;
    %  form clusters around present number of fault planes and compute
    %  lamda3 for each cluster instead of global variance J
    
    J=Copy_of_pcluster_using_lamda3(n0);
    
    max_lamda3 = max(lambda3);
 
    %fprintf('J = %g\n',J);
    %fprintf('Max lamda3 = %g\n',max_lamda3);
    
    
    %fprintf('Max lamda3 = %g\n',max_lamda3);
    
    if kj == 1
        max_lamda3OLD=max_lamda3;
    else
        max_lamda3NEW=max_lamda3;
        max_lamda3J=abs(max_lamda3OLD-max_lamda3NEW);
        max_lamda3OLD=max_lamda3NEW;
    end
    
    
%     fprintf('J= %g\n',J);
%     
%     if kj == 1;
%         JOLD=J;
%     else
%         JNEW=J;
%         DJ=abs(JOLD-JNEW);
%         JOLD=JNEW;
%     end
    
    %  Compute Cxy matrix, perform principal components analysis, create
    %  new fault planes for each cluster
    %  recalcfault(n0); % Seun moved it inside Copy_of_pcluster_using_lamda3(n0)
    
end

max_lamda3FINAL=max_lamda3NEW;

return


