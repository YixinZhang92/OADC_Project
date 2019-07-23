function read_catalog(infile,simul_tag)
%  read_catalog - read the hypocenters file to analyze

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc

fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

% plot the input data
% 
%figure;
Fig1 = figure('Name','Input hypocenters','Position', get(0, 'Screensize'));

plot3(xs,ys,zs,'o');
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

