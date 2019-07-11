function plot_strike_and_dip(grid_search_results)   

[m,~]=size(grid_search_results);
vp= 6.0; vs= 3.5; dens=2.7; z_displ=1; r_displ=1; t_displ=1;
line_color = 'r';
p_pol = [1.0 1.0 -1.0   0   -1.0 -1.0];

for i = 1:m
    best_strikes = grid_search_results(i,1);
    best_dips = grid_search_results(i,2);
    best_rakes = grid_search_results(i,3);

    [c_all,icontour_all, numstart_all, numfin_all] = ...
        calc_P_nodal_planes(best_strikes,best_dips,best_rakes,vp,vs,dens, z_displ, r_displ, t_displ);

    % plotting the P-wave nodal surface of the best focal mechanism
    % foceqarea_many_FM(c_all,icontour_all,numstart_all,...
    % numfin_all,p_pol, data_flag)
    % data_flag = 0 (do not plot polarity data)
    %           = 1 (otherwise)
    foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,...
        p_pol,1,line_color); hold on;
end