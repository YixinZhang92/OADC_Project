clear all; close all; clc;

data = textread('src-2.10');
xs_now = data(:,5);
ys_now = data(:,6);
zs_now = -data(:,7);

nhypos = length(xs_now);


analy = [xs_now ys_now zs_now zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1) ...
    zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1)  zeros(nhypos,1)]; 



strikes = 35;
dips = 60;

con=pi/180.;
strike=strikes.*con;
dip=dips.*con;

% figure
% scatter3(xs_now,ys_now,zs_now,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); %(Coordiante (0,0) is at the lower right corner)
% grid MINOR; hold on

% removing the mean for the rotation to be about origins. 
xs_now = xs_now - mean(xs_now); 
ys_now = ys_now - mean(ys_now);
zs_now = zs_now - mean(zs_now);

% figure
% scatter3(xs_now,ys_now,zs_now,'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% xlabel('Longitude (km)','FontSize',18)
% ylabel('Latitude (km)','FontSize',18)
% zlabel('depth (km)','FontSize',18)
% title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); %(Coordiante (0,0) is at the lower right corner)
% grid MINOR; hold on



% concatenating x,y and z
R=[xs_now' ; ys_now' ;zs_now'];



%************** rotate into strike direction **********************
%Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];

Rstrike=Dstrike*R;


rxp(1:nhypos) = Rstrike(1,1:nhypos);
ryp(1:nhypos) = Rstrike(2,1:nhypos);
rzp(1:nhypos) = Rstrike(3,1:nhypos);

%************** rotate into dip direction ********************
Rstrike(1,1:nhypos) = Rstrike(1,1:nhypos) - mean(Rstrike(1,1:nhypos));
Rstrike(2,1:nhypos) = Rstrike(2,1:nhypos) - mean(Rstrike(2,1:nhypos));
Rstrike(3,1:nhypos) = Rstrike(3,1:nhypos) - mean(Rstrike(3,1:nhypos));

% Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];

Rdip=Ddip*Rstrike;

% add bac the barycenter
Rdip(1,1:nhypos) = Rdip(1,1:nhypos) + mean(Rstrike(1,1:nhypos));
Rdip(2,1:nhypos) = Rdip(2,1:nhypos) + mean(Rstrike(2,1:nhypos));
Rdip(3,1:nhypos) = Rdip(3,1:nhypos) + mean(Rstrike(3,1:nhypos));

rxp(1:nhypos) = Rdip(1,1:nhypos);
ryp(1:nhypos) = Rdip(2,1:nhypos);
rzp(1:nhypos) = Rdip(3,1:nhypos);

% Take the depth to have a minimum of zero. Thi does not overright the
% original dataset.
rzp = rzp-min(rzp);



figure
scatter3(rxp,ryp,rzp,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
xlabel('Longitude (km)','FontSize',18)
ylabel('Latitude (km)','FontSize',18)
zlabel('depth (km)','FontSize',18)
title('Earthquake Distribution in the Charlevoix Seismic Zone','FontSize',18); %(Coordiante (0,0) is at the lower right corner)
grid MINOR; hold on



% Concatenate these result with the original dataet to kep track of each
% point. Sort the data based on the new depths.
analy(:,4) = rxp';
analy(:,5) = ryp';
analy(:,6) = rzp';


% Sorting in the new depth
analy = sortrows(analy,6);

% Stepwisely move through the data in depth and the number of EQs in each
% block. we need to specify the width of the block. 
st_array = 0:max(rzp+1); 




