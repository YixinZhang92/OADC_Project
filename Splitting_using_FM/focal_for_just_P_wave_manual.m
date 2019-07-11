clear all; close all; clc;

%% Inputs
strikes = 250;
dips = 60;
rakes = 165;

vp= 6; vs=3.5; dens= 2.7;

%p_pol = [p_amp  p_az        p_inc       weight_p_amp  p_index p_eps];
p_pol = [12930   271.0365    27.461801831175     5    1    -1;
         2704    232.12      55.0736556719227    5    2    -1;
         -3135   116.7252    60.4914228157975    5    3    -1;
         33.1    34.77872    60.7884642041669    5    4    -1;
         -1715   174.2638    75.762800373442     5    5    -1;
         -127.5  40.07827    80.1389504705361    5    6    -1;
         136.9   66.99249    82.416314722158     5    7    -1];

%%
z_displ=1; r_displ=1; t_displ=1;
[c_all,icontour_all, numstart_all, numfin_all] = ...
    calc_P_nodal_planes(strikes,dips,rakes,vp,vs,dens, z_displ, r_displ, t_displ);

% plotting the P-wave nodal surface of the best focal mechanism
% foceqarea_many_FM(c_all,icontour_all,numstart_all,...
% numfin_all,p_pol, data_flag)
% data_flag = 0 (do not plot polarity data)
%           = 1 (otherwise)
foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,...
    p_pol,1,'r'); hold on; shg