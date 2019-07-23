function [data,R] = create_syn_rnd_hypo_3D(nhypos,mu,sigma,outfile)
% clear all; close all; clc
% mu= 0;
% sigma=1;
% outfile='rnd_hypos_2D.txt';
% nhypos = 1000;

% Create some random data
R = normrnd(mu,sigma,[nhypos 1]);

phi = randn(nhypos,1)*2; %from 0 to 2pi
theta = randn(nhypos,1)*2; %from 0 to 2pi

x= R.*sin(theta).*cos(phi);
y= R.*sin(theta).*sin(phi);
z = R.*cos(theta);

data = [x y z];

%plot3(x,y,z,'ro'); shg

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:nhypos
    
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[x(kk) y(kk) z(kk)]);
    
end