function [rxp, ryp] = rand_hypos_2D(L,err_av,nhypos,strike,rt)

%  construct set of random hypocenters on a line with an additional
%  random error.  The line is specified by -L/2 < x < L/2 and
%  strike (degrees) and position.  nhypos is the number of hypocenters desired.
%  rt is a constant translation vector of the origin.

% %  outfile is the name of the output file of hypocenter locations
%  L=20; zerr_av=0.5; nhypos=200; strike= 70; rt=[0 -5];
%  outfile = 'rand_hypo_2D_1.txt';

con=pi/180.;
strike=strike.*con;

close all;

n=nhypos/4;
rx=[ (randperm(n) + randn(1,n))  (randperm(n) + randn(1,n)) ...
    -(randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n))] ;

% scale to L
rx=rx.*L./(2.0.*n);

% now construct an error in y
ry=randn(1,nhypos).*err_av;

% figure;
% plot(rx,ry,'o');
% axis equal

R=[rx ; ry];

% rotate into strike direction

Dstrike=[ sin(strike)  -cos(strike); cos(strike) sin(strike)];

Rstrike=Dstrike*R;

rxp(1:nhypos)=Rstrike(1,1:nhypos);
ryp(1:nhypos)=Rstrike(2,1:nhypos);

% figure;
% plot(rxp,ryp,'o');
% axis equal;
% title('rotation into strike');
% xlabel('X km');
% ylabel('Y km');

%  translate the final rotated collection of points
rxp=rxp+rt(1);
ryp=ryp+rt(2);

% figure;
% plot(rxp,ryp,'o');
% axis equal;
% title('Final translation of hypocenters');
% xlabel('X km');
% ylabel('Y km');
% grid on;

% % write synthetic data to outfile
% fid=fopen(outfile,'w');
% for kk=1:nhypos;
%     
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',[rxp(kk) ryp(kk) rzp(kk)]);
%     
% end;
