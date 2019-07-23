function [p_pol, sv_pol, sh_pol, vp, vs, dens] = real_polarities...
    (dist,az,depth, halfspace, vel_model, p_data, sv_data, sh_data,...
    path_to_SAC, path_to_travelt) 

% weight_P_amp,weight_S_amp, 
%
p_index = p_data(:,1); sv_index = sv_data(:,1); sh_index = sh_data(:,1);
p_amp = p_data(:,2); sv_amp = sv_data(:,2); sh_amp = sh_data(:,2);
weight_p_amp = p_data(:,3); weight_sv_amp = sv_data(:,3); weight_sh_amp = sh_data(:,3);

%   make a synthetic data set to test the focal grid search program
% We only extract dist and az where the P, SV or SH data are recorded.
[p_az, p_inc, p_eps, ~, ~, ~] = calc_inc_from_travelt...
    (dist(p_index),az(p_index),depth, vel_model,halfspace, 'p', path_to_SAC, path_to_travelt);         
         
% [sv_az, sv_inc, sv_eps, ~, ~, ~] = calc_inc_from_travelt...
%     (dist(sv_index),az(sv_index),depth, vel_model,halfspace, 'sv', path_to_SAC, path_to_travelt);
         
[s_az, s_inc, s_eps, vp, vs, dens] = calc_inc_from_travelt...
    (dist(sh_index),az(sh_index),depth, vel_model,halfspace, 'sh', path_to_SAC, path_to_travelt);         
% vp = 6.55;
% vs = 3.7;
% dens = 2.789;
%% In case the P, SV and SH are not recorded at all stations. 
% We removed stations with zero weight
p_pol  = [p_amp  p_az  p_inc  weight_p_amp  p_index p_eps];     %p_pol = p_pol(p_pol(:,4)~=0,:);
%sv_pol = [sv_amp sv_az sv_inc weight_sv_amp sv_index]; %sv_pol = sv_pol(sv_pol(:,4)~=0,:);
sv_pol = [sv_amp s_az s_inc weight_sv_amp sv_index s_eps]; %sv_pol = sv_pol(sv_pol(:,4)~=0,:);
sh_pol = [sh_amp s_az s_inc weight_sh_amp sh_index s_eps]; %sh_pol = sh_pol(sh_pol(:,4)~=0,:);