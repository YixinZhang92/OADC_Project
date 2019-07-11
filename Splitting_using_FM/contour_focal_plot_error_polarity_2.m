function contour_focal_plot_error_polarity_2(strikes, dips, rakes, fig_no, depth_now)
% clear all; close all; clc;
% 
% strikes = [80 80 78 78 78 78 76 106 100 80 80 80 102 216 80 76 78 78 80 78 30 120];
% dips = [30 30 30 30 30 30 30 32 28 30 30 30 26 76 28 36 30 30 32 30 22 32];
% rakes = [116 116 114 114 116 114 114 154 150 116 116 116 150 68 110 126 114 114 118 114 88 144];
% fig_no = 1

% N = length(strikes); az=linspace(0,360,1000); inc=linspace(0,90,300); P_contour = zeros(300,1000);
% T_contour = zeros(300,1000); I_contour = zeros(300,1000);

N = length(strikes); az=linspace(0,360,361); inc=linspace(0,90,91); P_contour = zeros(91,361);
T_contour = zeros(91,361); I_contour = zeros(91,361);


for i = 1:N
%   calculate the moment tensor for a point dislocation
m=dismom(strikes(i),dips(i),rakes(i));
%
%   calculate the P, T, and I vectors for this moment tensor
svec=stressvec(m);

ind = find(svec(:,1) < 0); svec(ind,1) = 360 + svec(ind,1);

P_coord(i,:) = round(svec(1,:),1);
I_coord(i,:) = round(svec(2,:),1);    
T_coord(i,:) = round(svec(3,:),1);    
  
v = abs(az - P_coord(i,1)); az_ind = find(v == min(v));
vi = abs(inc - P_coord(i,2)); inc_ind = find(vi == min(vi)); 
P_contour(inc_ind, az_ind) = P_contour(inc_ind, az_ind) + 1;

v = abs(az - T_coord(i,1)); az_ind = find(v == min(v));
vi = abs(inc - T_coord(i,2)); inc_ind = find(vi == min(vi)); 
T_contour(inc_ind, az_ind) = T_contour(inc_ind, az_ind) + 1;

v = abs(az - I_coord(i,1)); az_ind = find(v == min(v));
vi = abs(inc - I_coord(i,2)); inc_ind = find(vi == min(vi));
I_contour(inc_ind, az_ind) = I_contour(inc_ind, az_ind) + 1;

end

[x,y] = find(P_contour> 0); va = P_contour(P_contour> 0);
c= [az(y)' inc(x)' va];
cc_p = sortrows(c,3);

[x,y] = find(T_contour> 0); va = T_contour(T_contour> 0);
c= [az(y)' inc(x)' va];
cc_t = sortrows(c,3);

[x,y] = find(I_contour> 0); va = I_contour(I_contour> 0);
c= [az(y)' inc(x)' va];
cc_i = sortrows(c,3);

PTIfocsphere(cc_p, cc_i, cc_t, fig_no, depth_now)

end


function PTIfocsphere(cc_p, cc_i, cc_t, fig_no, depth_now)
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
fig = figure(fig_no);
%figure('color','w','name','P Wave Nodal Surfaces');
hold on;
daspect([1 1 1]);
axis off;
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

[m,~]=size(cc_p);
color_line = [fliplr(linspace(0.1,0.8,max(cc_p(:,3))-1)) 0]; %1- (1: max(cc_p(:,3)))/max(cc_p(:,3));

for i = 1:m
    svec = [cc_p(i,1:2)];
    
    %   Plot P stress axes
    %
    rp=sqrt(2.0)*radius*sin(svec(1,2)*con/2.);
    yct=rp*cos(svec(1,1)*con);
    xct=rp*sin(svec(1,1)*con);
    mark='o'; %fcol='r';
    %
    val_p = color_line(cc_p(i,3));
    plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[1,val_p,val_p]);  
end
%

[m,~]=size(cc_i);
color_line = [fliplr(linspace(0.1,0.8,max(cc_i(:,3))-1)) 0]; %1- (1: max(cc_p(:,3)))/max(cc_p(:,3));

for i = 1:m
    svec = [cc_i(i,1:2)];
    
    %   Plot I stress axes
    %
    rp=sqrt(2.0)*radius*sin(svec(1,2)*con/2.);
    yct=rp*cos(svec(1,1)*con);
    xct=rp*sin(svec(1,1)*con);
    mark='o'; %fcol='r';
    %
    val_i = color_line(cc_i(i,3));
    plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[val_i,1,val_i]);  
end


[m,~]=size(cc_t);
color_line = [fliplr(linspace(0.1,0.8,max(cc_t(:,3))-1)) 0]; %1- (1: max(cc_p(:,3)))/max(cc_p(:,3));

for i = 1:m
    svec = [cc_t(i,1:2)];
    
    %   Plot T stress axes
    %
    rp=sqrt(2.0)*radius*sin(svec(1,2)*con/2.);
    yct=rp*cos(svec(1,1)*con);
    xct=rp*sin(svec(1,1)*con);
    mark='o'; %fcol='r';
    %
    val_t = color_line(cc_t(i,3));
    plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[val_t,val_t,1]);  
end

%     %
%     %   Plot I stress axes
%     %
%     rp=sqrt(2.0)*radius*sin(svec(2,2)*con/2.);
%     yct=rp*cos(svec(2,1)*con);
%     xct=rp*sin(svec(2,1)*con);
%     mark='o';
%     val_i = color_line(cc_i(i,3));
%     plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[val_i,1,val_i]); 
%     
%     %
%     %   Plot T stress axes
%     %
%     rp=sqrt(2.0)*radius*sin(svec(3,2)*con/2.);
%     yct=rp*cos(svec(3,1)*con);
%     xct=rp*sin(svec(3,1)*con);
%     mark='o'; 
%     val_t = color_line(cc_t(i,3));
%     plot(xct,yct,'marker',mark,'markeredgecolor','k','markerfacecolor',[val_t,val_t,1]); 
%     
    
    
    
    
    
    
    
%colorbar
markersize = 8;
textfontsize = 14;
%
title(['Error due to Polarity (depth = ' num2str(depth_now) ' km)'], 'FontSize', 20);
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

    
print(fig, ['fig_err_' num2str(fig_no) '.png'], '-dpng') 
   
   
end

