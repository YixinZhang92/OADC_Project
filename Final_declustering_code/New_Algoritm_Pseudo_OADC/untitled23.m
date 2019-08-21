clear all; close all; clc;

rng(5,'twister');
X = mvnrnd([0 0 0], [1 .2 .7; .2 1 0; .7 0 1],50);
plot3(X(:,1),X(:,2),X(:,3),'bo');
grid on;
maxlim = max(abs(X(:)))*1.1;
axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
axis square
view(-9,12);


[coeff,score,roots] = pca(X);
basis = coeff(:,1:2);
normal = coeff(:,3);


pctExplained = roots' ./ sum(roots);

[n,p] = size(X);
meanX = mean(X,1);
Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';
residuals = X - Xfit;


error = abs((X - repmat(meanX,n,1))*normal);
sse = sum(error.^2);


[xgrid,ygrid] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),2), ...
                         linspace(min(X(:,2)),max(X(:,2)),2));
zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
%h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);

% xgrid = xgrid';
% ygrid = ygrid';
% zgrid = zgrid';


xv=[xgrid(1,1) xgrid(2,1) xgrid(2,2) xgrid(1,2)]; 
yv=[ygrid(1,1) ygrid(2,1) ygrid(2,2) ygrid(1,2)]; 
zv=[zgrid(1,1) zgrid(2,1) zgrid(2,2) zgrid(1,2)]; 

hold on;
for k=1%:n0
    fill3(xv(1:4),yv(1:4),zv(1:4),'w','FaceAlpha',0.6,'FaceColor','k');
end






hold on
above = (X-repmat(meanX,n,1))*normal < 0;
below = ~above;
nabove = sum(above);
X1 = [X(above,1) Xfit(above,1) nan*ones(nabove,1)];
X2 = [X(above,2) Xfit(above,2) nan*ones(nabove,1)];
X3 = [X(above,3) Xfit(above,3) nan*ones(nabove,1)];
plot3(X1',X2',X3','-', X(above,1),X(above,2),X(above,3),'o', 'Color',[0 .7 0]);
nbelow = sum(below);
X1 = [X(below,1) Xfit(below,1) nan*ones(nbelow,1)];
X2 = [X(below,2) Xfit(below,2) nan*ones(nbelow,1)];
X3 = [X(below,3) Xfit(below,3) nan*ones(nbelow,1)];
plot3(X1',X2',X3','-', X(below,1),X(below,2),X(below,3),'o', 'Color',[1 0 0]);

% hold off
% maxlim = max(abs(X(:)))*1.1;
% axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
% axis square
% view(-9,12);


