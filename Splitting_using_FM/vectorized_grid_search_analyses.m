function [grid_search_results] = vectorized_grid_search_analyses...
    (strike_incr,dip_incr,rake_incr,p_pol,sv_pol,sh_pol,vp,vs,dens,...
    w_pol, w_ratio, use_sv_amp, use_sv_pol, z_displ, r_displ, t_displ)
% This function performs a grid search method to determine the best focal
% mechanism that fits the P, SV and SH wave polarities.
%
% Inputs: strike_incr,dip_incr,rake_incr = strike, dip and rake increments
%         p_pol,sv_pol,sh_pol = P-, SV- and SH-wave polarity information
%         vp,vs,dens = P- and S-wave velocities, and density
%
% Outputs: grid_search_results = strike, dip, rake, Percent_acc
%          mm = number of correct(best) focal mechanisms based on P,SV and SH
%          polarities
%
% strike_incr = 5;
% dip_incr = 5;
% rake_incr = 5;
% vp = 6;
% vs = 3.5;
% dens = 2.7;

max_strike = (360-strike_incr);
max_dip = 90;
max_rake = (360-rake_incr);

strike = 0:strike_incr:max_strike;
dip = 0:dip_incr:max_dip;
rake = 0:rake_incr:max_rake;     
% 
% fprintf('version 3\n'); tic
% determining moment tensor components
[m1,m2,m3,m4,m5,m6]=vect_dismom_v3(strike,dip,rake);

% calculating wave amplitudes based on moment tensor components and source param
% Each strike, dip and rake combination is contained in each row of prad,
% svrad and shrad and each column is for each station.
% 
% In order words, the dimension of prad, svrad and shrad is (Nai) x (Nd*Ns*Nr) where 
% Nai - number of stations       Nd - number of dip andgles
% Ns - Number of strike angles   Nr - Number of rake angles
prad = vect_prad_vect(m1,m2,m3,m4,m5,m6,p_pol,vp,vs,dens, z_displ);
[svrad, svrad_r] = vect_svrad_vect(m1,m2,m3,m4,m5,m6,sv_pol,vp,vs,dens, r_displ);
shrad = vect_shrad_vect(m1,m2,m3,m4,m5,m6,sh_pol,vp,vs,dens, t_displ);
%
%% GRIDSEARCH USING P, SV AND SH POLARITY DATA
% converting amplitude data to pola and removing NaN from the polarity data
% I am using only P and SH polarities for real data. They are more reliable.
p_pol_data = convert_data_to_polarity(p_pol(:,1)); 
p_pol_data(isnan(p_pol_data))=0;
sv_pol_data = convert_data_to_polarity(sv_pol(:,1)); 
sv_pol_data(isnan(sv_pol_data))=0;
sh_pol_data = convert_data_to_polarity(sh_pol(:,1)); 
sh_pol_data(isnan(sh_pol_data))=0;
%
% Applying the weight of the amplitudes
% I want to first normalize the weights for it to range from 0 to 1.
weight_p_amp = p_pol(:,4);      w_p_norm = weight_p_amp/5;
weight_sv_amp = sv_pol(:,4);    w_sv_norm = weight_sv_amp/5;
weight_sh_amp = sh_pol(:,4);    w_sh_norm = weight_sh_amp/5;

w_p_pol_data = p_pol_data.*w_p_norm;
w_sv_pol_data = sv_pol_data.*w_sv_norm;
w_sh_pol_data = sh_pol_data.*w_sh_norm;

% performing dot product of observed and calculated polarities
prad_pol = convert_data_to_polarity(prad); prad_pol(isnan(prad_pol))=0;
svrad_pol = convert_data_to_polarity(svrad); svrad_pol(isnan(svrad_pol))=0;
shrad_pol = convert_data_to_polarity(shrad); shrad_pol(isnan(shrad_pol))=0;

% use_sv_pol = 0 -- do not use sv polarity
% use_sv_pol = 1 -- use sv polarity
if use_sv_pol == 0
    p_dot_prod = prad_pol' * w_p_pol_data;
    sh_dot_prod = shrad_pol' * w_sh_pol_data;

    C = (p_dot_prod+ sh_dot_prod);
    N = sum(w_p_norm) + sum(w_sh_norm);
    % This was added so that low weight do not further decrease the accuracy.
    %N = length(p_pol(:,1))+length(sh_pol(:,1));
    
elseif use_sv_pol == 1
    p_dot_prod = prad_pol' * w_p_pol_data;
    sv_dot_prod = svrad_pol' * w_sv_pol_data;
    sh_dot_prod = shrad_pol' * w_sh_pol_data;

    C = (p_dot_prod+ sv_dot_prod + sh_dot_prod);
    N = sum(w_p_norm) + sum(w_sh_norm) + sum(w_sh_norm); 
    % This was added so that low weight do not further decrease the accuracy.
    %N = length(p_pol(:,1))+length(sh_pol(:,1))+length(sv_pol(:,1));

end

ns=length(strike); nr=length(rake); nd=length(dip);
%
% Calculating percentag of accuracy. If the synthetic and real data are the
% same, the dot product would be N. So, the percent_acc = 100%
Percent_acc_pol = (C./N)*100; 
%save Percentpol1.mat Percent_acc_pol;
%
%% GRIDSEARCH USING SV/P, SH/P AND SV/SH AMPLITUDE RAIOS DATA
% calculating atan of amplitude ratios and sum them for every strike, dip
% and rake combinations.

% use_sv_amp = 0 -- do not use sv at all
% use_sv_amp = 1 -- use only |sv|/sh
% use_sv_amp = 2 -- use only |sv|/p
% use_sv_amp = 3 -- use both |sv|/sh and |sv|/p

% In the amplitude ratio gridsearch, we need only the radial component of 
% svrad since that is what is measured by the user. However, we use an absolute
% value of svrad to plot the S-wave polarization vector of the focalsphere.
svrad = svrad_r;

if use_sv_amp == 1
    % I want to use the absolute value of SV amplitude. This is because of the
    % the polarity changes that could occur due to wave conversions.
    svrad = abs(svrad); sv_pol(:,1) = abs(sv_pol(:,1));
    %svrad = (svrad); sv_pol(:,1) = (sv_pol(:,1));
    % 
    Percent_acc_sh_p = calc_percent_acc_amp_ratios(shrad,prad,sh_pol,p_pol);
    Percent_acc_sv_sh = calc_percent_acc_amp_ratios(shrad,svrad,sh_pol,sv_pol);
    %
    Percent_acc_amp_ratios = (Percent_acc_sh_p' + Percent_acc_sv_sh')/2;

elseif use_sv_amp == 2
    svrad = abs(svrad); sv_pol(:,1) = abs(sv_pol(:,1));
    %svrad = (svrad); sv_pol(:,1) = (sv_pol(:,1));
    % 
    Percent_acc_sv_p = calc_percent_acc_amp_ratios(svrad,prad,sv_pol,p_pol);
    Percent_acc_sh_p = calc_percent_acc_amp_ratios(shrad,prad,sh_pol,p_pol);
    %
    Percent_acc_amp_ratios = (Percent_acc_sh_p' + Percent_acc_sv_p')/2;
    
elseif use_sv_amp == 3
    svrad = abs(svrad); sv_pol(:,1) = abs(sv_pol(:,1));
    %svrad = (svrad); sv_pol(:,1) = (sv_pol(:,1));
    % 
    Percent_acc_sv_p = calc_percent_acc_amp_ratios(svrad,prad,sv_pol,p_pol);
    Percent_acc_sh_p = calc_percent_acc_amp_ratios(shrad,prad,sh_pol,p_pol);
    Percent_acc_sv_sh = calc_percent_acc_amp_ratios(shrad,svrad,sh_pol,sv_pol);
    %
    Percent_acc_amp_ratios = (Percent_acc_sh_p' + Percent_acc_sv_p' +...
    Percent_acc_sv_sh')/3;
    
elseif use_sv_amp == 0
    % 
    Percent_acc_sh_p = calc_percent_acc_amp_ratios(shrad,prad,sh_pol,p_pol);
    %
    Percent_acc_amp_ratios = Percent_acc_sh_p';
end
 

%%
Percent_acc = ((w_pol*Percent_acc_pol)) + (w_ratio*Percent_acc_amp_ratios);
%
%save Percent1.mat Percent_acc;

Percent_acc = reshape(Percent_acc,[nd ns nr]);

% determine the unique values in the percent of accuracy of the focal
% mechanisms. We need only the best 2 mechanisms.
a = unique(Percent_acc);
% a
% vp
% vs
% dens
[dip_r0,strike_r0,rake_r0] = ind2sub(size(Percent_acc),find(Percent_acc == a(end)));
[dip_r1,strike_r1,rake_r1] = ind2sub(size(Percent_acc),find(Percent_acc == a(end-1)));

% Creating Output table
Nstr0 = length(strike_r0); Nstr1 = length(strike_r1);

grid_search_results(:,1)= [strike(strike_r0)';strike(strike_r1)']; % strike
grid_search_results(:,2)= [dip(dip_r0)';dip(dip_r1)']; % dip
grid_search_results(:,3)= [rake(rake_r0)';rake(rake_r1)']; % rake
grid_search_results(:,4)= [a(end)*ones(Nstr0,1);a(end-1)*ones(Nstr1,1)]; % percent_acc
