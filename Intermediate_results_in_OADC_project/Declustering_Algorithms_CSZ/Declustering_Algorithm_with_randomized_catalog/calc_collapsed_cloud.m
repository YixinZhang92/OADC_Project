function [x_colap, y_colap, z_colap] = calc_collapsed_cloud(x, y, z, uncertainty, n_iter)

% x=xs; y= ys; z = zs;
% uncertainty = 1;
%  n_iter = 5;

N = length(x);
x_temp(1:N)=0.0;
y_temp(1:N)=0.0;
z_temp(1:N)=0.0;


for kk = 1: n_iter
    for k=1:N
     
        x0=x(k);
        y0=y(k);
        z0=z(k);

        % calculate distance from reference source location
        dist=sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2);

        % determine the hypocenters within 1 km of the object hypocenter
        x_c = x(dist<=uncertainty);
        y_c = y(dist<=uncertainty);
        z_c = z(dist<=uncertainty);

        % calc centroid of the determined hypocenters
%         x_temp(k) = mean(x_c);
%         y_temp(k) = mean(y_c);
%         z_temp(k) = mean(z_c);        
        
        x_temp(k) = x0 + 0.61803*(mean(x_c) - x0);
        y_temp(k) = y0 + 0.61803*(mean(y_c) - y0);
        z_temp(k) = z0 + 0.61803*(mean(z_c) - z0);
        
%         (x0 - mean(x_c))
    end
    
    x = x_temp;
    y = y_temp;
    z = z_temp; 
    
end

x_colap = x;
y_colap = y;
z_colap = z; 

