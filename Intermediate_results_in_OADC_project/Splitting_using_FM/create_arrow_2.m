function create_arrow_2(x,y,length,angle,arrowcolor)
%  [p0_bar, p1_bar, ah_1_bar,ah_2_bar] = 
% clear all; clc; close all;
% x= 0; y = 0; length = 5; angle = 30; arrowcolor = 'r';
%
% This function creates a vertical arrow, rotate clockwisely by angle and
% translate it to point (x,y).
%
% creating an horizontal arrow with unit length, enlarge it based on
% given length.
half_len = length/2;
p0=-[0 half_len];
p1=[0 half_len];
ah_1 = p1 + 0.5*half_len*[sind(30) -cosd(30)];
ah_2 = p1 + 0.5*half_len*[-sind(30) -cosd(30)];
%
% Applying rotating and translation. The rotation is clockwise.
rot_mat = [cosd(angle) sind(angle); -sind(angle) cosd(angle)];
%
p0_bar = (rot_mat*p0') + [x; y];
p1_bar = (rot_mat*p1') + [x; y];
ah_1_bar = (rot_mat*ah_1') + [x; y];
ah_2_bar = (rot_mat*ah_2') + [x; y];
%
% how to project the angles on the lower hemisphere
plot([p0_bar(1) p1_bar(1)],[p0_bar(2) p1_bar(2)],arrowcolor); hold on;
plot([p1_bar(1) ah_1_bar(1)],[p1_bar(2) ah_1_bar(2)],arrowcolor); hold on;
plot([p1_bar(1) ah_2_bar(1)],[p1_bar(2) ah_2_bar(2)],arrowcolor); hold on; 