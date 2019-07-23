function [x_colap, y_colap, z_colap, dist_to_ref_loc] = calc_collapsed_cloud(x, y, z, uncertainty, n_iter, skim)

% x=xs; y= ys; z = zs;
% skim = 'JS_';
% uncertainty=1; n_iter =2;

N = length(x);
x_temp(1:N)=0.0;
y_temp(1:N)=0.0;
z_temp(1:N)=0.0;

x_ini = x;
y_ini = y;
z_ini = z;

dist_to_ref_loc = zeros(n_iter,N);

x_colap = zeros(n_iter,N);
y_colap = zeros(n_iter,N);
z_colap = zeros(n_iter,N);


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
        if strcmp(skim,'JS')
            x_movi = 0.61803*(mean(x_c) - x0);
            y_movi = 0.61803*(mean(y_c) - y0);
            z_movi = 0.61803*(mean(z_c) - z0);
            
        elseif strcmp(skim,'Nichol_etal')
            dist_c = dist(dist<=uncertainty);
            weight_c = (1/sqrt(2*pi))*exp(-dist_c.^2/2);

            x_movi = 0.61803*(sum(weight_c.*x_c)/sum(weight_c) - x0);
            y_movi = 0.61803*(sum(weight_c.*y_c)/sum(weight_c) - y0);
            z_movi = 0.61803*(sum(weight_c.*z_c)/sum(weight_c) - z0);
            
        else
            disp('***************************')
            disp('Enter the rght skim to use')
            disp('***************************')
            break;
        end
        
        x_temp(k) = x0 + x_movi;
        y_temp(k) = y0 + y_movi;
        z_temp(k) = z0 + z_movi;
        
    end
    
    x = x_temp;
    y = y_temp;
    z = z_temp; 
  
    % calculate distance of new location from reference/initial location
    dist_to_ref_loc(kk,:)=sqrt((x-x_ini).^2 + (y-y_ini).^2 + (z-z_ini).^2);
    
    x_colap(kk,:) = x;
    y_colap(kk,:) = y;
    z_colap(kk,:) = z; 
    
end
