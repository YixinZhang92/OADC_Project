function datplot_with_colors(xs,ys,zs,n0,xv,yv,zv,picname,simul_tag)

% Seun editted datplot.m to included simul_tag, and save figures with the tag
% as the filename. He also included figure title and adjusted fontsize.


% plot a rendition of the data and the planes that fit the data
%figure('Name',picname);
Fig1 = figure('Name',picname,'Position', get(0, 'Screensize'));

hold on;
plot3(xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor','k');


colors=[0. 1 1; 
    0.5 0.5 0.8;
    0.7   0   0;
    1   0.8 0.8;
    0.5 0.8 0.5;
    0.5 0.8 0.5;
    0.7 0.9 0.5];

for k=1:n0;
    fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.6,'FaceColor',colors(k,:));
end;
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




legend('Eqs','Fault 1','Fault 2','Fault 3','Fault 4','Fault 5','Fault 6',...
    'Fault 7','Location','southoutside','Orientation','horizontal');




%xlim([10 70]);

% Printing figure to file
fig_filename = [simul_tag '.' picname '.png'];
F1    = getframe(Fig1);
imwrite(F1.cdata, fig_filename, 'png')
savefig(Fig1,[fig_filename(1:end-4) '.fig'])

return;