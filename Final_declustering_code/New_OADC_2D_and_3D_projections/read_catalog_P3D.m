function read_catalog_P3D(infile,simul_tag,opt)
%  read_catalog - read the hypocenters file to analyze

global xc yc zc vec_plane xb_old yb_old zb_old orig_xs orig_ys orig_zs N Nc

fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
orig_xs(1:N)=data(1:N,1);
orig_ys(1:N)=data(1:N,2);
orig_zs(1:N)=data(1:N,3);

if opt == 1
    % plot the input data
    % 
    %figure;
    Fig1 = figure('Name','Input hypocenters','Position', get(0, 'Screensize'));

    plot3(orig_xs,orig_ys,orig_zs,'o');
    axis equal;
    grid on;
    title('Input hypocenters');
    xlabel('X km');
    ylabel('Y km');
    zlabel('Z km');

    set(gca, 'fontsize', 18);

    % Printing figure to file
    fig_filename = [simul_tag '.input.hypos.png'];
    F1    = getframe(Fig1);
    imwrite(F1.cdata, fig_filename, 'png')
    savefig(Fig1,[fig_filename(1:end-4) '.fig'])
end
end

