function Percent_acc_sv_p = calc_percent_acc_amp_ratios(svrad,prad,sv_pol,p_pol)
%% GRIDSEARCH BY AMPLITUDE RAIOS DATA
% calculating atan of amplitude ratios and sum them for every strike, dip
% and rake combinations.
% sv and p here is just an example. This function can work for other
% amplitude ratios.
% Note: sv_p_calc has the same dimensions as the prad. i.e. (Nai) x (Nd*Ns*Nr)
%
%  WE USE ONLY STATIONS WITH BOTH P AND SV AMPLITUDES
%condition_p_sv = (p_pol(:,5) == sv_pol(:,5)); % station index must the same
condition_p = find(ismember(p_pol(:,5),sv_pol(:,5),'rows'));
condition_sv = find(ismember(sv_pol(:,5),p_pol(:,5),'rows'));

%% calculated angles 
% Syntax: an = atan_S_wave_pol(data_x,data_y) % station index must the same
svrad_now = svrad(condition_sv,:); 
prad_now = prad(condition_p,:);
%
sv_p_calc = atan_S_wave_pol(svrad_now,prad_now); %save sv_p_calc.mat sv_p_calc
%
%size(sv_p_calc)
%% observed angles
p_pol_amp = p_pol(condition_p,1);
sv_pol_amp = sv_pol(condition_sv,1);
%


% % The original size of the angles is (Nai) x 1 but the size of 
% % calculated angles is (Nai) x (Nd*Ns*Nr). So, we need to repeat the column
% % array (Nd*Ns*Nr) times for us to do matrix subraction and compare each
% % column of the calculated with observed angle (column).



%  siz_data = size(sv_p_calc); no_of_FM = siz_data(2);
% sv_p_obs = repmat(atan_S_wave_pol(sv_pol_amp,p_pol_amp),1,no_of_FM); % observed % save sv_p_obs.mat sv_p_obs
% %
% %size(sv_pol_amp)
% %phi_sv_p = sum(abs(sv_p_calc - sv_p_obs),1); %save phi_sv_p.mat phi_sv_p
% %phi_sv_p(isnan(phi_sv_p))=0;
% %
% % determining the smallest angle difference.
% phi_sv_p1 = abs(sv_p_calc - sv_p_obs); %save phi_sv_p.mat phi_sv_p

sv_p_obs = atan_S_wave_pol(sv_pol_amp,p_pol_amp);
%bsxfun(@minus, matrix, column)
phi_sv_p1 = abs(bsxfun(@minus, sv_p_calc, sv_p_obs)); 

phi_sv_p1(isnan(phi_sv_p1))=0;
indices = find(abs(phi_sv_p1)>=180);
phi_sv_p1(indices) = 360 - phi_sv_p1(indices);
% save phi_sv_p1.mat phi_sv_p1

% Before we sum the angle differences, we want to apply the weight at each
% station. We want to apply the minimum weigths of the wave at each
% station and add all the angle differences for each focal mechanism.

weight_p = p_pol(condition_p,4)/5; 
weight_sv = sv_pol(condition_p,4)/5;
weight_com = [weight_p weight_sv];

weight_min = min(weight_com,[],2) ;
%size(weight_min)
phi_sv_p =  weight_min' * phi_sv_p1;
%phi_sv_p = sum(phi_sv_p1,1);c
%
% maximum deviation at each station
siz_data = size(sv_p_calc);


phi_sum_sv_p = 179*sum(weight_min); 
% This was added so that low weight do not further decrease the accuracy.
%phi_sum_sv_p = 179*siz_data(1); 


% phi_sum_sv_p = 359*siz_data(1); 
%
% percentage of accuracy
Percent_acc_sv_p = (1 - (phi_sv_p/phi_sum_sv_p))*100;