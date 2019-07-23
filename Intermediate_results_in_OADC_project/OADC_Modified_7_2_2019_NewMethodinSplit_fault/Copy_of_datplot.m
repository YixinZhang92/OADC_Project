function Copy_of_datplot(xs,ys,zs,n0,xv,yv,zv,picname,simul_tag)

% plot a rendition of the data and the planes that fit the data
Fig1 = figure('Name',picname,'Position', get(0, 'Screensize'));

hold on;
plot3(xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
for k=1:n0
    fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end
hold off;
axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
grid on;
view(3);
title(picname);
set(gca, 'fontsize', 18);

% Printing figure to file
fig_filename = [simul_tag '.' picname '.png'];
F1    = getframe(Fig1);
imwrite(F1.cdata, fig_filename, 'png')
savefig(Fig1,[fig_filename(1:end-4) '.fig'])

return;