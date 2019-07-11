function focal_error_analysis_focaldepth(strike_incr,dip_incr,rake_incr,...
    p_pol,sv_pol,sh_pol, w_pol, w_ratio, dist, az, depth, halfspace, ...
    vel_model, depth_error, use_sv_amp, use_sv_pol, path_to_SAC, ...
    path_to_travelt, z_displ, r_displ, t_displ)
%
% This function determines the error analysis for focal.
% It determines the effect of focal depth on the calculated focal
% mechanisms. At each focal depth, it loops over each observed amplitudes and changes its polarities.
% In each loop, this function performs grid search to determine sets of
% focal mechanisms that fit the observed polarities (P, SV and SH
% polarities) and amplitude ratios (SV/P, SH/P and SV/SH). It determines
% the modal focal mechanism at each focal depth and plot the strike, dip
% and rake as a function of focal depth. It also plots histogram of strike,
% dip and rake at the given focal depth as we change the amplitude polarities.
% 
% To test this function, uncomment the section and comment the function
%
% clear all; close all; clc;
% strike_incr=5;dip_incr=5;rake_incr=5;
% [p_pol, sv_pol, sh_pol] = load_correct_data;
% w_pol = 0.7; w_ratio = 0.3;
% vel_model =[6.08   3.51   2.732    0.0    1; ...
%             6.25   3.60   2.756    6.0   1; 
%             6.55   3.70   2.789    12.0   1;
%             7.2    4.20   2.90    40.0   1;
%             8.0    4.60   3.2    40.0001   2;
%             8.0    4.60   3.2    200.0   2]; 
%         
% depth = 24.5; depth_error = 2; halfspace = 0;       
% dist = [3.5265 20 34.6410 9.3262 16.7820 28.5630 5.3590 11.5470 ...
%         14.0042 20 9.3262 20 34.6410 7.2794];
% az = [10 15. 25. 50. 80. 105. 158. 190. 210. 250. 270. 300. 320. 333.];
         
% ------------ Error analysis begins -------------------
[l_p, ~] = size(p_pol); [l_sv, ~] = size(sv_pol); 
pol = [p_pol(:,1); sv_pol(:,1); sh_pol(:,1)]; iter = length(pol);

weight_p_amp = p_pol(:,4); weight_sv_amp = sv_pol(:,4); weight_sh_amp = sh_pol(:,4);
p_index = p_pol(:,5); sv_index = sv_pol(:,5); sh_index = sh_pol(:,5);

depth_min = depth - depth_error;
if depth_min < 0; depth_min = 0.; end

depth_max = depth + depth_error; jj = 0;
depth_range = depth_min:depth_max;

% determine updated azimuth and incidence angle at each focal depth
% and flip the polarities of each amplitude at each iteration.
for j = 1:length(depth_range)
    jj = jj + 1; depth_now = depth_range(jj);
    
    fprintf('\n')
    fprintf(['Performing error analysis at depth = ', num2str(depth_now) ' km \n'])
    
    % determine azimuth and inc angle using travelt program
    [p_az, p_inc,p_eps, ~, ~, ~] = calc_inc_from_travelt...
        (dist(p_index),az(p_index),depth_now, vel_model,halfspace,...
        'p', path_to_SAC, path_to_travelt);         
         
%     [sv_az, sv_inc,eps, ~, ~, ~] = calc_inc_from_travelt...
%         (dist(sv_index),az(sv_index),depth_now, vel_model,halfspace,...
%         'sv', path_to_SAC, path_to_travelt);
         
    [s_az, s_inc, s_eps, vp, vs, dens] = calc_inc_from_travelt...
        (dist(sh_index),az(sh_index),depth_now, vel_model,halfspace,...
        'sh', path_to_SAC, path_to_travelt);         
  
    % concantenate for p_pol, sv_pol and sh_pol
    p_pol_depth = [p_pol(:,1) p_az p_inc weight_p_amp p_index p_eps];
    sv_pol_depth = [sv_pol(:,1) s_az s_inc weight_sv_amp sv_index s_eps];
    sh_pol_depth = [sh_pol(:,1) s_az s_inc weight_sh_amp sh_index s_eps];
    
    for i = 1:(iter+1)
        % store the updated azimuths and inc angles for the error analysis.
        % Note that the amplitude do not change throughout the error analysis.
        p_pol_err = p_pol_depth; sv_pol_err = sv_pol_depth; sh_pol_err = sh_pol_depth;
        pol_error = pol; % initializing the original amplitudes

        % i=1 corresponds to the original dataset
        if i > 1; pol_error(i-1) = pol_error(i-1)*(-1); end 
       
        % split the amplitudes to the corresponding waves
        p_pol_err(:,1) = pol_error(1:l_p);
        sv_pol_err(:,1) = pol_error(l_p+1:l_p+l_sv);
        sh_pol_err(:,1) = pol_error(l_p+l_sv+1:end);
        
        % perform grid search (Vectorized version)
        grid_search_results = vectorized_grid_search_analyses...
            (strike_incr,dip_incr,rake_incr,p_pol_err,sv_pol_err,...
            sh_pol_err, vp,vs,dens, w_pol, w_ratio, use_sv_amp, ...
            use_sv_pol, z_displ, r_displ, t_displ);
  
        % extract strike, dip and rake from the results from grid search analysis
        strikes(i) = grid_search_results(1,1);
        dips(i) = grid_search_results(1,2); 
        rakes(i) = grid_search_results(1,3);
        bestFM_acc(i) = grid_search_results(1,4);

    end

    % determines the model of strikes, dips and rakes at each focal depth
    strike_mode(jj) = mode(strikes);
    dip_mode(jj) = mode(dips);
    rake_mode(jj) = mode(rakes);
   
%     save strikes.mat strikes
%     save dips.mat dips 
%     save rakes.mat rakes 
    
    
    fig_no = j + 100;
    contour_focal_plot_error_polarity_2(strikes, dips, rakes, fig_no, depth_now);
    
    
%     ESPP = 1e-10;
%     if ((depth - depth_now) <= ESPP)
%         figure (4)
%         subplot(1,3,1)
%         hist(strikes,(min(strikes):strike_incr:max(strikes)));
%         title(['Histogram of Strike Values (Strike = ', num2str(mode(strikes)),' )'])
% 
%         subplot(1,3,2)
%         hist(dips,(min(dips):dip_incr:max(dips))); 
%         title(['Histogram of Dip Values (Dip = ', num2str(mode(dips)),' )'])
% 
%         subplot(1,3,3)
%         hist(rakes,(min(rakes):rake_incr:max(rakes))); 
%         title(['Histogram of Rake Values (Rake = ', num2str(mode(rakes)),' )'])
%     end
   
end   

% Plots
fig_err_w_depth = figure (j+101);

subplot(3,1,1)
plot(depth_range,strike_mode,'r'); grid MINOR;
title('Effect of Focal Depth on Strike')

subplot(3,1,2)
plot(depth_range,dip_mode,'b'); grid MINOR;
title('Effect of Focal Depth on Dip')

subplot(3,1,3)
plot(depth_range,rake_mode,'g'); grid MINOR;
title('Effect of Focal Depth on Rake')

axis off;
print(fig_err_w_depth, 'fig_err_w_depth.png', '-dpng') 
   
shg    

% --------- END OF CODE ----------------






















% %% 
% figure (1)
% subplot(1,3,1)
% %plot(strike_1,bestFM_acc,'ro');
% hist(strikes,(min(strikes):strike_incr:max(strikes)));
% title(['Histogram of Strike Values (Strike = ', num2str(mode(strikes)),' )'])
% 
% subplot(1,3,2)
% %plot(dip_1,bestFM_acc,'bo');
% hist(dips,(min(dips):dip_incr:max(dips))); 
% title(['Histogram of Dip Values (Dip = ', num2str(mode(dips)),' )'])
% 
% subplot(1,3,3)
% %plot(rake_1,bestFM_acc,'go');
% hist(rakes,(min(rakes):rake_incr:max(rakes))); 
% title(['Histogram of Rake Values (Rake = ', num2str(mode(rakes)),' )'])
% 
% figure (2)
% subplot(1,3,1)
% plot(strikes,bestFM_acc,'ro');
% title(['Histogram of Strike Values (Strike = ', num2str(mode(strikes)),' )'])
% 
% subplot(1,3,2)
% plot(dips,bestFM_acc,'bo');
% title(['Histogram of Dip Values (Dip = ', num2str(mode(dips)),' )'])
% 
% subplot(1,3,3)
% plot(rakes,bestFM_acc,'go');
% title(['Histogram of Rake Values (Rake = ', num2str(mode(rakes)),' )'])