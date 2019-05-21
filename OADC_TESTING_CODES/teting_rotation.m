clear all; close all; clc;
global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc

L=10; W=5; zerr_av=0.1; nhypos=200;strike=30;dip=90;rake=0;rt=[0 0 -5];
outfile = 'rand_hypo_seun.txt';

rand_hypos(outfile,L,W,zerr_av,nhypos,strike,dip,rake,rt)

read_catalog(outfile);



% Random hypocenters have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike

% rotate into rake direction

nhypos = length(xs);
strike = 30;%pos is anticlockwise
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






%Drake=[ cos(rake)  -sin(rake) 0 ; sin(rake) cos(rake) 0 ; 0 0 1];
Drake=[ cos(rake) 0 -sin(rake); 0 1 0; sin(rake) 0 cos(rake)];

Rrake=Drake*Rdip;

rxp(1:nhypos)=Rrake(1,1:nhypos);
ryp(1:nhypos)=Rrake(2,1:nhypos);
rzp(1:nhypos)=Rrake(3,1:nhypos);

figure;
plot3(rxp,ryp,rzp,'o');
axis equal;
title('Rotation into rake');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');



