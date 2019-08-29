function datplot_2D(xs,ys,n0,xv,yv,picname,simul_tag)

% Seun editted datplot.m to included simul_tag, and save figures with the tag
% as the filename. He also included figure title and adjusted fontsize.


% plot a rendition of the data and the planes that fit the data
%figure('Name',picname);
Fig1 = figure('Name',picname,'Position', get(0, 'Screensize'));

hold on;
plot(xs,ys,'o','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
for k=1:n0
    %fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
    plot(xv(k,1:2),yv(k,1:2),'r');
end

hold off;
axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
grid on;

title(picname);
set(gca, 'fontsize', 18);

% Printing figure to file
fig_filename = [simul_tag '.' picname '.png'];
F1    = getframe(Fig1);
imwrite(F1.cdata, fig_filename, 'png')
savefig(Fig1,[fig_filename(1:end-4) '.fig'])

return;