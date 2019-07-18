

clear all; close all; clc;


xst =0; yst=0; zst=0;

Lt = 10;
Wt = 2;


xvt = [(-Wt/2) (Wt/2) (Wt/2) (-Wt/2)];
yvt = [(-Lt/2) (-Lt/2) (Lt/2) (Lt/2)];
zvt = [0 0 0 0];

kmint = 1; simul_tagt = 'simul';

strike = 35; dip = 65; rake = 0;

%%
con=pi/180.;
strike=strike.*con;
dip=dip.*con;
rake=rake.*con;

% Fault corner coordinates have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike

% rotate into rake direction

R=[xvt ; yvt ;zvt];
[~,n] = size(R);

Drake=[ cos(rake)  -sin(rake) 0 ; sin(rake) cos(rake) 0 ; 0 0 1];

Rrake=Drake*R;

xvt(1:n)=Rrake(1,1:n);
yvt(1:n)=Rrake(2,1:n);
zvt(1:n)=Rrake(3,1:n);

% picnamet = 'rotate.into.rake';
% 
% datplot(xst,yst,zst,kmint,xvt,yvt,zvt,picnamet,simul_tagt);


% rotate into dip direction

%Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
Ddip=[ cos(dip) 0 sin(dip); 0 1  0 ; -sin(dip) 0 cos(dip)];

Rdip=Ddip*Rrake;

xvt(1:n)=Rdip(1,1:n);
yvt(1:n)=Rdip(2,1:n);
zvt(1:n)=Rdip(3,1:n);

% picnamet = 'rotate.into.dip';
% 
% datplot(xst,yst,zst,kmint,xvt,yvt,zvt,picnamet,simul_tagt);


% rotate into strike direction

%Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
Dstrike=[ cos(strike)  sin(strike) 0 ; -sin(strike) cos(strike) 0 ; 0 0 1];

Rstrike=Dstrike*Rdip;

xvt(1:n)=Rstrike(1,1:n);
yvt(1:n)=Rstrike(2,1:n);
zvt(1:n)=Rstrike(3,1:n);

picnamet = 'rotate.into.strike';

datplot(xst,yst,zst,kmint,xvt,yvt,zvt,picnamet,simul_tagt);













% 
% 
% figure;
% plot3(rxp,ryp,rzp,'o');
% axis equal;
% title('Rotation into dip');
% xlabel('X km');
% ylabel('Y km');
% zlabel('Z km');
% 
% 
% % rotate into strike direction
% 
% Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
% 
% Rstrike=Dstrike*Rdip
% 
% rxp(1:nhypos)=Rstrike(1,1:nhypos);
% ryp(1:nhypos)=Rstrike(2,1:nhypos);
% rzp(1:nhypos)=Rstrike(3,1:nhypos);
% 
% figure;
% plot3(rxp,ryp,rzp,'o');
% axis equal;
% title('rotation into strike');
% xlabel('X km');
% ylabel('Y km');
% zlabel('Z km');
% 
% %  translate the final rotated collection of points
% rxp=rxp+rt(1);
% ryp=ryp+rt(2);
% rzp=rzp+rt(3);
% 
% figure;
% plot3(rxp,ryp,rzp,'o');
% axis equal;
% title('Final translation of hypocenters');
% xlabel('X km');
% ylabel('Y km');
% zlabel('Z km');
% grid on;
% 
% % write synthetic data to outfile
% fid=fopen(outfile,'w');
% for kk=1:nhypos;
%     
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',[rxp(kk) ryp(kk) rzp(kk)]);
%     
% end;
% 
% 
% %%
% 
% 
% 
% 
% 
% 
% 
% datplot(xst,yst,zst,kmint,xvt,yvt,zvt,picnamet,simul_tagt);
