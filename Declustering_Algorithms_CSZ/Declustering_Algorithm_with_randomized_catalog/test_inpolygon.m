clear all; close all; clc;

L = linspace(0,2.*pi,6);
xv = cos(L)';
yv = sin(L)';

rng default
xq = randn(250,1);
yq = randn(250,1);

[in,on] = inpolygon(xq,yq,xv,yv);

numel(xq(in))

numel(xq(~in))


figure

plot(xv,yv) % polygon
axis equal

hold on
plot(xq(in),yq(in),'r+') % points inside
plot(xq(~in),yq(~in),'bo') % points outside
hold off