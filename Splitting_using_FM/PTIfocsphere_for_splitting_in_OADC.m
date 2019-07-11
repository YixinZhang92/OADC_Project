function PTIfocsphere_for_splitting_in_OADC(grid_search_results)
%
%   produce a figure of the PTI for each focal mechanism determined 
% from the grid search analyses on an equal area, lower hemisphere, projection
%
%   Input:  grid_search_results = strike, dip, rake, Percent_acc
%*********************************************************
%
%   create the outside circle of the projection
%
con=pi/180.;
radius=1.0;
a=linspace(0,360,1000)*con;
yc=radius*cos(a);
xc=radius*sin(a);
%
%figure;
%figure('color','w','name','P Wave Nodal Surfaces');
hold on;
daspect([1 1 1]);
%axis off;
plot(xc,yc,'linewidth',1.5,'color',[0.5,0.5,0.5]);
%
%   plot N,E,S,W and center tick marks
%
for it=1:4
    az=90*(it-1)*con;
    r1=radius*0.97;
    r2=radius*1.03;
    xt=radius*1.1*sin(az);
    yt=radius*1.1*cos(az);
    plot([r1*sin(az) r2*sin(az)],[r1*cos(az) r2*cos(az)]);
    if it == 1
        text(xt,yt,'N','horizontalalignment','center', 'FontSize',14);
    end
    if it == 2
        text(xt,yt,'E','horizontalalignment','center', 'FontSize',14);
    end
    if it == 3
        text(xt,yt,'S','horizontalalignment','center', 'FontSize',14);
    end
    if it == 4
        text(xt,yt,'W','horizontalalignment','center', 'FontSize',14);
    end
end
plot([-0.03 0.03],[0 0]);
plot([0 0],[-0.03 0.03]);
%
% determine the unique values in the percent correctness of the focal
% mechanisms. We need only the best 3 mechanisms.
%a = unique(grid_search_results(:,4));
%
[m,~]=size(grid_search_results);
for i = m:-1:1
    strike = grid_search_results(i,1);
    dip = grid_search_results(i,2);
    rake = grid_search_results(i,3);
    %quality = grid_search_results(i,4);
    %
    %   calculate the moment tensor for a point dislocation
    m=dismom(strike,dip,rake);
    %
    %   calculate the P, T, and I vectors for this moment tensor
    %
    svec=stressvec(m);
    %
    %   Plot P stress axes
    %
    rp=sqrt(2.0)*radius*sin(svec(1,2)*con/2.);
    yct=rp*cos(svec(1,1)*con);
    xct=rp*sin(svec(1,1)*con);
    mark='o'; %fcol='r';
    %
%     if perc >= a(end)
         plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[1,0,0]); 
%     elseif perc >= a(end-1)
%         plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[1,0.8,0.8]);
%     end
    %
    %   Plot I stress axes
    %
 
    rp=sqrt(2.0)*radius*sin(svec(2,2)*con/2.);
    yct=rp*cos(svec(2,1)*con);
    xct=rp*sin(svec(2,1)*con);
    mark='o'; 
%     if perc >= a(end)
         plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[0,1,0]); 
%     elseif perc >= a(end-1)
%         plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[0.8,1,0.8]);
%     end
%     %
    %   Plot T stress axes
    %
    rp=sqrt(2.0)*radius*sin(svec(3,2)*con/2.);
    yct=rp*cos(svec(3,1)*con);
    xct=rp*sin(svec(3,1)*con);
    mark='o'; 
%     if perc >= a(end)
         plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[0,0,1]);   
%     elseif perc >= a(end-1)
%         plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[0.8,0.8,1]);
%     end   
end
%
markersize = 8;
textfontsize = 14;
%
title('PTI stress axes', 'FontSize', 20);
%  Put polarity info onto plot
plot(0.9,1.2,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor','r');
text(1.0,1.2,'P-axis', 'FontSize', textfontsize);
plot(0.9,1.05,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor','b');
text(1.0,1.05,'T-axis', 'FontSize', textfontsize);
plot(0.9,0.90,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor','g');
text(1.0,0.90,'I-axis', 'FontSize', textfontsize);
%
format long g
%
% plot(0.9,-1.05,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor',[1,0,0]);
% text(1.0,-1.05,strcat(num2str(round(a(end)*100)/100),'%'), 'FontSize', textfontsize);
% plot(0.9,-0.9,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor',[1,0.8,0.8]);
% text(1.0,-0.9,strcat(num2str(round(a(end-1)*100)/100),'%'), 'FontSize', textfontsize);
% %
hold off; 
return;