X=[xs' ys' zs'];

% compute the covariance matrix for this cluster
Cxy=cov(X,0);
[V,D]=eig(Cxy);
basis1=V(1:3,2)
basis2=V(1:3,3)
vec_plane=V(1:3,1);






[coeff,score,roots] = pca(X);
basis = coeff(:,1:2)
normal = coeff(:,3);


