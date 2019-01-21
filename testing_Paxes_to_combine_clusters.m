clear all; close all; clc;

index = [186;70;77;214;213;193;172;179;226;17;187;233;248];

no = [124;79;53;13;12;5;3;3;3;2;1;1;1];

strike = [130;45;50;150;150;135;120;125;160;10;130;165;175];

dip = [45;90;90;45;30;45;45;45;15;30;60;15;30];

rake = 0;

for i = 1: length(index)

%   calculate the moment tensor for a point dislocation
m=dismom(strike(i),dip(i),rake);
%
%   calculate the P, T, and I vectors for this moment tensor
%
svec=stressvec(m);
%

az_inc(i,1:5) = [no(i) strike(i) dip(i) svec(1,1) svec(1,2)];

end

az_inc(az_inc(:,4)<0,4) = az_inc(az_inc(:,4)<0,4) + 360;

% %   Plot P stress axes
% %
% rp=sqrt(2.0)*radius*sin(svec(1,2)*con/2.);
% yct=rp*cos(svec(1,1)*con);
% xct=rp*sin(svec(1,1)*con);
% mark='o'; %fcol='r';