function [x_colap, y_colap, nSD_moved] = calc_collapsed_cloud_2D...
    (x, y, nSD, sigma_x, sigma_y, n_iter, skim)

x=xs; y= ys; %z = zs;
% skim = 'JS_';
% uncertainty=1; n_iter =2;

N = length(x);
x_temp(1:N)=0.0;
y_temp(1:N)=0.0;
%z_temp(1:N)=0.0;

x_ini = x;
y_ini = y;
%z_ini = z;

nSD_moved = zeros(n_iter,N);

x_colap = zeros(n_iter,N);
y_colap = zeros(n_iter,N);
%z_colap = zeros(n_iter,N);

% determine the ellipse scale
CI = erf(nSD/sqrt(2));
s = chi2inv(CI,nSD);

for kk = 1%: n_iter
    for k=1:N
        
        x0=x_ini(k);
        y0=y_ini(k);
        %z0=z(k);

        % calculate distance from reference source location
        s_dist=(((x-x0)/sigma_x).^2 + ((y-y0)/sigma_y).^2);% + ((z-z0)/sigma_z).^2);

        % determine the hypocenters within 1 km of the object hypocenter
        x_c = x(s_dist<=s);
        y_c = y(s_dist<=s);
        %z_c = z(s_dist<=s);
   
        
        % calc centroid of the determined hypocenters
        if strcmp(skim,'JS')
            x_movi = 0.61803*(mean(x_c) - x0);
            y_movi = 0.61803*(mean(y_c) - y0);
            %z_movi = 0.61803*(mean(z_c) - z0);
            
        elseif strcmp(skim,'Nichol_etal')
            % calculate distance from reference source location
            dist=sqrt((x-x0).^2 + (y-y0).^2);% + (z-z0).^2);
        
            dist_c = dist(s_dist<=s);
            weight_c = (1/sqrt(2*pi))*exp(-dist_c.^2/2);

            x_movi = 0.61803*(sum(weight_c.*x_c)/sum(weight_c) - x0);
            y_movi = 0.61803*(sum(weight_c.*y_c)/sum(weight_c) - y0);
            %z_movi = 0.61803*(sum(weight_c.*z_c)/sum(weight_c) - z0);
            
        else
            disp('***************************')
            disp('Enter the rght skim to use')
            disp('***************************')
            break;
        end
        
        x_temp(k) = x0 + x_movi;
        y_temp(k) = y0 + y_movi;
        %z_temp(k) = z0 + z_movi;
        
    end
    
    x = x_temp;
    y = y_temp;
    %z = z_temp; 
  
%     % calculate distance of new location from reference/initial location
%     dist_to_ref_loc(kk,:)=sqrt((x-x_ini).^2 + (y-y_ini).^2 + (z-z_ini).^2);
    
    % calculate distance from reference source location in terms of s
    sini_dist=(((x-x_ini)/sigma_x).^2 + ((y-y_ini)/sigma_y).^2);% + ((z-z0)/sigma_z).^2);

    % convert s to number of SD moved from theobject hypocenter
    CI_ini = chi2cdf(sini_dist, nSD);
    nSD_moved = sqrt(2)* erfinv(CI_ini);
        
    x_colap(kk,:) = x;
    y_colap(kk,:) = y;
    %z_colap(kk,:) = z; 
    
end




        
        
        
        
        
        
% chisquare_val = sqrt(s);%2.4477;
% theta_grid = linspace(0,2*pi);
% phi = 0;%angle;
% % X0=avg(1);
% % Y0=avg(2);
% a=chisquare_val*sqrt(1);
% b=chisquare_val*sqrt(1);
% 
% % the ellipse in x and y coordinates 
% ellipse_x_r  = a*cos( theta_grid );
% ellipse_y_r  = b*sin( theta_grid );
% 
% %Define a rotation matrix
% R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
% 
% %let's rotate the ellipse to some angle phi
% r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
% 
% % Draw the error ellipse
% figure; plot(x,y,'ro'); hold on; plot(x_c,y_c,'bo'); hold on;
% plot(r_ellipse(:,1)+ x0,r_ellipse(:,2) + y0,'-')
% hold on;
% plot(x0,y0,'go')




