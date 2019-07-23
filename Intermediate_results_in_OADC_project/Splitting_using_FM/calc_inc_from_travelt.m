function [az_now, inc_now, eps, vp, vs, dens] = calc_inc_from_travelt...
    (x_intp,az,depth, vel_model,halfspace, wave, path_to_SAC, path_to_travelt)
%
vp_vel = vel_model(:,1); vs_vel = vel_model(:,2);
dens_vel = vel_model(:,3); depth_vel = vel_model(:,4);

vp = interp1(depth_vel,vp_vel,depth,'linear');
vs = interp1(depth_vel,vs_vel,depth,'linear');
dens = interp1(depth_vel,dens_vel,depth,'linear');
%
% create in file for travelt program
%create_in_file_for_travelt(wave, x_intp, depth);
create_in_file_for_travelt_halfspace(wave, x_intp, depth, vel_model, halfspace);
%
% Initializing travelt and SAC
setenv('SACAUX',path_to_SAC);
path = path_to_travelt; 

%
prefix = ['CSZ_travelt_',wave];
%
% Running travelt program
%! /Users/oluwaseunfadugba/Documents/Waveform_Modeling_Now/time/travelt ... 
% < CSV_travelt.in > CSV_travelt.out
eval(sprintf('%s%s%s%s%s%s%s','!',path,'/travelt <',prefix,'.in > ',prefix,'.out'));
!rm -r data; mkdir data; mv *CSZ_L* data; rm -f *~
%
% Reading the results from travelt program
fnames_x_t = dir('data/*x-t*'); %CSZ_Layered_sh.x-t.h001.001
fnames_x_p = dir('data/*x-p*');
numfids = length(fnames_x_t); % number of vtk files
% 
%fnames_x_t.name
index = 0;
for K = 1:numfids % iterate over the files
   filename_x_t = (fnames_x_t(K).name); % extract the filename
   filename_x_p = (fnames_x_p(K).name); % extract the filename
   
   % I want to remove SP wave. It arrives earlier than SV or SH waves.
   if strcmp(wave,'sv') || strcmp(wave,'sh')
       if strcmp(filename_x_t(end-4:end),'1.001')...
           || strcmp(filename_x_t(end-4:end),'2.001') || ...
           strcmp(filename_x_t(end-4:end),'3.001') || ...
           strcmp(filename_x_t(end-4:end),'4.001') 
        
           continue;
       end
       %wave
       index = index + 1;
   else
       index = index + 1;
   end 
   
   [x,t]=readsac(sprintf('%s%s','data/',filename_x_t));
   [~,p]=readsac(sprintf('%s%s','data/',filename_x_p));
   %
   filename_x_t_now{index} = filename_x_t;
   %filename_x_t
   % Interpolate the functions
   % I removed extrapolation because extrapolation of headwaves gives an 
   % incorrect traveltime at smaller distances   
   t_intp(index,:) = interp1(x,t,x_intp,'linear');  % index was previously K
   p_intp(index,:) = interp1(x,p,x_intp,'linear'); 
  
%    %% Plot x_t and x_p 
%    subplot(2,1,1)
%    grid MINOR;
%    plot(database(N,x_ind),database(N,t_ind)); grid MINOR; hold on; 
%     title([wave,' wave Travel Time Vs Distance'])
%    
%    subplot(2,1,2)
%    grid MINOR;
%    plot(database(N,x_ind),database(N,p_ind)); hold on; grid MINOR;
%    title([wave,' wave Ray Parameter Vs Distance'])
 end
%  xlim([0 300]);
%  legend(fnames_x_t.name); 
%
 %t_intp
 %filename_x_t_now{4}
 
 %t_intp
 %p_intp
 
for n_x = 1:length(x_intp)
    t_now(n_x) = min(t_intp(:,n_x));
    % index of the first arrival
    index_of_min = find(t_intp(:,n_x) == min(t_intp(:,n_x)));
    % extract the ray parameter at the first arrival using the index
    p_now(n_x,1) = p_intp(index_of_min, n_x); 
    % extract the filename for first arrival
    %up_or_downgoing{n_x,1} = fnames_x_t(index_of_min).name; This is wrong
    %because it still takes the index_of_min from the list of all the
    %original files. It should be the list where S-P conversion have been
    %removed.
    
    up_or_downgoing{n_x,1} = filename_x_t_now(index_of_min);
    %
    %fnames_x_t(index_of_min).name
    % determine if the ray is upgoing
    %if strfind(fnames_x_t(index_of_min).name,'.001')  This is also wrong
    %the same reason above.
    
    temp_f = (filename_x_t_now(index_of_min));
    temp_file = temp_f{1};
    
    if strcmp(temp_file(end-3:end),'.001')
        azimuth_plus(n_x,1) = 0; % formerly 180. But I don' want to add 180 
        % deg here. I will add it only for plotting on the focal sphere.
        eps(n_x,1) = -1;
    else
        azimuth_plus(n_x,1) = 0;
        eps(n_x,1) = 1;
    end
end
%
%p_now
if strcmp(wave,'p'); 
    inc_ang = asin(p_now*vp)*180*7/22;
else
    inc_ang = asin(p_now*vs)*180*7/22;
end
%
if halfspace ==1
    azimuth_plus = zeros(length(az),1);
end
%
azimuth_now = az + azimuth_plus';
%
az_now = azimuth_now';
inc_now = inc_ang;
%
% Remove temporary files
!rm -r data; rm -f CSZ_travelt_*

%eps
%azimuth_plus




%clear all; close all; clc;
% %% Inputs
% x_intp = [3.5265 20 34.6410 9.3262 16.7820 28.5630 5.3590 11.5470 ...
%     14.0042 20 9.3262 20 34.6410 7.2794];
% az = [10 15. 25. 50. 80. 105. 158. 190. 210. 250. 270. 300. 320. 333.];
% amplitude = [1  1.  1.  1.  1.  1.   1.   1.   1.   1.   1.   1.   1.   1.];
% 
% wave = 'sv';
% depth = 20;
% vel_model = [6.08    3.51   2.737    0.0; ... % Halfspace model
%              6.08    3.51   2.737    40.0]; 
%
% setenv('SACAUX','/usr/local/sac/bin');
% path = '/Users/oluwaseunfadugba/Documents/Waveform_Modeling_Now/time';