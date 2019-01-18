


% Initialize the array of x,y,z,norm,index
% 
% for every strike and dip, extract points that falls within 4km, and check if they are more than 3
%     
% calc eig and check if lamda3 < 2
% 
% determine strike and dip, and check if they are les than 10 from the plan been tested
% 
% determine W and check if > 5
% 
% calc normalized factor and assign it to all points. If normalized factor exisit, choose the smaller with it index


close all; clear all; clc;
infile = 'testdata.txt';
global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale


%***************** Read Catalog of Hypocenters ****************************
read_catalog(infile);

% Random hypocenters have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike

% rotate into rake direction
tic
nhypos = length(xs);
strike = 45;%pos is anticlockwise
dip = 90;
rake = 0;

con=pi/180.;
strike=strike.*con;
dip=dip.*con;
rake=rake.*con;


xs = xs - mean(xs);
ys = ys - mean(ys);
zs = zs - mean(zs);

R=[xs ; ys ;zs];


% rotate into strike direction

%Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];

Rstrike=Dstrike*R;

rxp(1:nhypos)=Rstrike(1,1:nhypos);
ryp(1:nhypos)=Rstrike(2,1:nhypos);
rzp(1:nhypos)=Rstrike(3,1:nhypos);

figure;
plot3(rxp,ryp,rzp,'o');
axis equal;
title('rotation into strike');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

% rotate into dip direction
Rstrike(1,1:nhypos) = Rstrike(1,1:nhypos)-mean(Rstrike(1,1:nhypos));
Rstrike(2,1:nhypos) = Rstrike(2,1:nhypos)-mean(Rstrike(2,1:nhypos));
Rstrike(3,1:nhypos) = Rstrike(3,1:nhypos)-mean(Rstrike(3,1:nhypos));

% Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
%Ddip=[ cos(dip) 0 sin(dip); 0 1  0 ; -sin(dip) 0 cos(dip)];
Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];

Rdip=Ddip*Rstrike;

Rdip(1,1:nhypos) = Rdip(1,1:nhypos)+mean(Rstrike(1,1:nhypos));
Rdip(2,1:nhypos) = Rdip(2,1:nhypos)+mean(Rstrike(2,1:nhypos));
Rdip(3,1:nhypos) = Rdip(3,1:nhypos)+mean(Rstrike(3,1:nhypos));

rxp(1:nhypos)=Rdip(1,1:nhypos);
ryp(1:nhypos)=Rdip(2,1:nhypos);
rzp(1:nhypos)=Rdip(3,1:nhypos);

figure;
plot3(rxp,ryp,rzp,'o');
axis equal;
title('Rotation into dip');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

rzp = rzp-min(rzp);

analy = [xs; ys; zs; rxp; ryp; rzp]';
analy=sortrows(analy,6);

st_array= 0:max(rzp+1) ; width= 1;
for i = 1: length(st_array)
    st = st_array(i);
    aa(i) = length(analy(analy(:,6)>=st & analy(:,6)<=(st+width),6));
end

aa(aa<=20) = 0;

plot(st_array, aa); shg

toc