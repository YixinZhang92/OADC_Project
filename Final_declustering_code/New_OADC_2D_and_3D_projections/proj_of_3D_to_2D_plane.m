function [pxs,pys] = proj_of_3D_to_2D_plane(xs,ys,zs,az,el)

X=[xs' ys' zs'];

U = [cosd(az), -sind(az), 0];
V = [-sind(az)*sind(el), -cosd(az)*sind(el), -cosd(el)];
n = [sind(az)*cosd(el), cosd(az)*cosd(el), -sind(el)];

R = [U' V' n'];
projX=R\X';
pxs = projX(1,:);
pys = -projX(2,:);

% Fig1 = figure('Name','Projected hypocenters','Position', get(0, 'Screensize'));
% subplot(1,2,1)
% plot3(xs,ys,zs,'o');axis equal; view(-az, el);grid MINOR
% 
% subplot(1,2,2)
% plot(pxs,pys,'o');grid MINOR; shg

% 
% outfile = [simul_tag 'proj.' num2str(az) '.' num2str(el) '.txt'];
% 
% % write synthetic data to outfile
% fid=fopen(outfile,'w');
% for kk=1:N
%     
%     fprintf(fid,'%12.5f %12.5f\n',[pxs(kk) pys(kk)]);
%     
% end

