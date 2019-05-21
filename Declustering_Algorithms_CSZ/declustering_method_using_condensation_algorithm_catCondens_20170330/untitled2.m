xb = mean(xs)
yb = mean(ys)
zb = mean(zs)



 dist=sqrt((xs-xb).^2 + (ys-yb).^2 + (zs-zb).^2);
 
 
 [cdf_rnd,vr]=ecdf(dist);

 figure;
 plot(vr(2:end),cdf_rnd-cdf_rnd(1),'-k','LineWidth',2);hold on;
% semilogx(vr(2:end),diff(cdf_rnd),'-k','LineWidth',2);hold on;
%plot([log(V05) log(V05)],[0 1],'-r');
%semilogx([V05 V05],[0 1],'-r','LineWidth',2);
hold off;

title('Cumulative Distribution Functions');
xlabel('volume (km^3)');
ylabel('Probability');
legend('Random','Catalog');
grid MINOR;
set(gca, 'fontsize', 18);

figure
ecdfhit(dist)

figure
[counts, binCenters] = hist(dist, 256); % Use 256 bins.
plot(binCenters, counts);



y2 = pdf('Chisquare',3)

y2 = pdf('Normal',dist,0,4)
