function plot_focal_figures...
    (grid_search_results, panel,p_pol,sv_pol,sh_pol, vp, vs, dens, ...
    z_displ, r_displ, t_displ)
%
% extracting strike, dip and rakes, and the fitting accuracy for the best 
% focal mechanisms from the grid search method.
a = unique(grid_search_results(:,4));
condition(:,1) = (grid_search_results(:,4) == a(end));
condition(:,2) = (grid_search_results(:,4) == a(end-1));

line_color=['r','b'];
if strcmp(panel,'p_panel')
    % Plot the PTI of the FMs on a focal sphere
    PTIfocsphere(grid_search_results);
    
elseif strcmp(panel,'sv_panel')   
    for ii = [2 1];
        best_strikes = grid_search_results(condition(:,ii),1);
        best_dips = grid_search_results(condition(:,ii),2);
        best_rakes = grid_search_results(condition(:,ii),3);

        [c_all,icontour_all, numstart_all, numfin_all] = ...
            calc_P_nodal_planes(best_strikes,best_dips,best_rakes,vp,vs,dens, z_displ, r_displ, t_displ);

        % plotting the P-wave nodal surface of the best focal mechanism
        % foceqarea_many_FM(c_all,icontour_all,numstart_all,...
        % numfin_all,p_pol, data_flag)
        % data_flag = 0 (do not plot polarity data)
        %           = 1 (otherwise)
        foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,...
            p_pol,1,line_color(ii)); hold on;
    end
    
elseif strcmp(panel,'sh_panel')
    for ii = [2 1];
        best_strikes = grid_search_results(condition(:,ii),1);
        best_dips = grid_search_results(condition(:,ii),2);
        best_rakes = grid_search_results(condition(:,ii),3);

        [n_bst_str, ~] = size(best_strikes);
        for j = 1:n_bst_str;
            % calculated angles 
            % fprintf('version 3\n'); tic
            % determining moment tensor components
            %best_strikes,best_dips,best_rakes
            [m1,m2,m3,m4,m5,m6]=vect_dismom_v3...
                (best_strikes(j),best_dips(j),best_rakes(j));

            % calculating wave amplitudes based on moment tensor components 
            % and source param. Each strike, dip and rake combination is 
            % contained in each row of prad, svrad and shrad and each column 
            % is for each station.
            % 
            % In order words, the dimension of prad, svrad and shrad is 
            % (Nai) x (Nd*Ns*Nr) where 
            % Nai - number of stations         Nd - number of dip andgles
            % Ns - Number of strike angles     Nr - Number of rake angles          
            svrad = vect_svrad_vect(m1,m2,m3,m4,m5,m6,sv_pol,vp,vs,dens, r_displ);
            shrad = vect_shrad_vect(m1,m2,m3,m4,m5,m6,sh_pol,vp,vs,dens, t_displ);

            
%             % This function updates the azimuth at each station based on the direction
%             % of wave propagation. 
%             % az = az + 180 for upgoing waves
%             % az = az       for downgoing waves
%             sv_pol = update_az(sv_pol);
%             sh_pol = update_az(sh_pol);
%             
            

            %  WE USE ONLY STATIONS WITH BOTH SV AND SH AMPLITUDES
            % azimuth must the same
            condition_sv_sh = (sv_pol(:,2) == sh_pol(:,2)); 
            sv_pol_amp = sv_pol(condition_sv_sh,1);
            sh_pol_amp = sh_pol(condition_sv_sh,1);
            az = sv_pol(condition_sv_sh,2);
            
            % update_az do not work because of the for loop keeps increasing the az.
            s_eps = sv_pol(condition_sv_sh,6); 
            I = (s_eps == -1); az(I) = az(I)+180;
            
            inc = sv_pol(condition_sv_sh,3);

            % Syntax: an = atan_S_wave_pol(data_x,data_y)
            % Because of the geometry of SV and SH, we have to switch the x and
            % y axes. It has to be x-SH and y-SV which is equivalent to SH/SV
            % rather than SV/SH ratio. switched. save sv_p_calc.mat sv_p_calc
            rot_angle_cal_1 = az + atan_S_wave_pol(shrad,svrad);
            %
            %focus = [az inc atan_S_wave_pol(shrad,svrad) shrad svrad]; %save focus.mat focus
            %
            % cal S-wave pol. with maximum accuracy
            plot_Swave_pola_foceqarea(svrad,shrad, az, inc, ...
                rot_angle_cal_1,line_color(ii)) 
        end
    end
    %
    % The original size of the angles is (Nai) x 1 but the size of 
    % calculated angles is (Nai) x (Nd*Ns*Nr). So, we need to repeat the column
    % array (Nd*Ns*Nr) times for us to do matrix subraction and compare each
    % column of the calculated with observed angle (column).
    %
    % rot_angle_obs_now = repmat(atan_S_wave_pol(sv_pol_amp,sh_pol_amp),1,...
    % no_of_FM); % observed % save sv_p_obs.mat sv_p_obs
    % Because of the geometry of SV and SH, we have to switch the x and
    % y axes. It has to be x-SH and y-SV which is equivalent to SH/SV
    % rather than SV/SH ratio. % switched. observed % save sv_p_obs.mat sv_p_obs
    rot_angle_obs_now = az + atan_S_wave_pol(sh_pol_amp,sv_pol_amp); 
  
    % superimpose the S-wave polarization on the focal sphere observed S-wave pol.
    plot_Swave_pola_foceqarea(sv_pol_amp,sh_pol_amp, az, inc, rot_angle_obs_now,'k')

    markersize = 8; textfontsize = 14;

    title('S Wave Polarization*', 'FontSize', 20);
    %  Put polarity info onto plot
    plot(0.9,1.2,'MarkerSize',markersize,'marker','o','markeredgecolor',...
        'k','markerfacecolor','k');
    text(1.0,1.2,'Observed Spol', 'FontSize', textfontsize);
    plot(0.9,1.05,'MarkerSize',markersize,'marker','o','markeredgecolor',...
        'k','markerfacecolor','r');
    text(1.0,1.05,'Calc. Spol_1', 'FontSize', textfontsize);
    plot(0.9,0.90,'MarkerSize',markersize,'marker','o','markeredgecolor',...
        'k','markerfacecolor','b');
    text(1.0,0.90,'Calc. Spol_2', 'FontSize', textfontsize);
    hold on;
    
    for ii = [2 1];
        best_strikes = grid_search_results(condition(:,ii),1);
        best_dips = grid_search_results(condition(:,ii),2);
        best_rakes = grid_search_results(condition(:,ii),3);

        [c_all,icontour_all, numstart_all, numfin_all] = ...
            calc_P_nodal_planes(best_strikes,best_dips,best_rakes,vp,vs,...
                                dens, z_displ, r_displ, t_displ);

        % plotting the P-wave nodal surface of the best focal mechanism
        % foceqarea_many_FM(c_all,icontour_all,numstart_all,...
        % numfin_all,p_pol, data_flag)
        % data_flag = 0 (do not plot polarity data)
        %           = 1 (otherwise)
        foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,...
            p_pol,0,line_color(ii)); hold on
    end
    title('S Wave Polarization*', 'FontSize', 20);
end