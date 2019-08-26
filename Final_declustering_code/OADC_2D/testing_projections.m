clear all; close all; clc;
global xs ys zs N

infile = 'Simul.1_hypos.txt';
simul_tag = 'Simul.5Faults.OADC';
read_catalog(infile,simul_tag);

figure
for az=40%:5:360
for el =25%0:5:360
    
theta = az;
phi = el;

% U = [-sind(theta), cosd(theta), 0];
% V = [cosd(theta)*sind(phi), sind(theta)*sind(phi), cosd(phi)];
% Center = [cosd(theta)*cosd(phi), sind(theta)*cosd(phi), sind(phi)];

U = [cosd(theta), -sind(theta), 0];
V = [-sind(theta)*sind(phi), -cosd(theta)*sind(phi), -cosd(phi)];
Center = [sind(theta)*cosd(phi), cosd(theta)*cosd(phi), -sind(phi)];

X = [xs' ys' zs'];

for i=1:N
x(i) = dot(X(i,:) - Center, U);
y(i) = dot(X(i,:) - Center, V);

word(i,:) = Center + U * x(i) + V * y(i);

end

plot(x,y,'ro'); %xlim([-15 25]);ylim([-45 0]); shg

% subplot(1,2,1)
% plot(x,y,'ro'); %xlim([-15 25]);ylim([-45 0]); shg
% subplot(1,2,2)
% plot3(word(:,1),word(:,2),word(:,3),'ro'); shg
pause(0.5)
end
end

outfile = ['proj.' num2str(az) '.' num2str(el) '.txt'];

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:N
    
    fprintf(fid,'%12.5f %12.5f\n',[x(kk) y(kk)]);
    
end







% 
% V(1,1) = cosd(el)*cosd(az);
% V(2,1) = cosd(el)*sind(az);
% V(3,1) = sind(el);
% 
% X = [xs' ys' zs'];
% 
% X_2d = X - (X*V)*V';
% 
% figure
% subplot(1,2,1)
% plot3(X(:,1), X(:,2), X(:,3),'bo'); view(az,el)
% subplot(1,2,2)
% plot3(X_2d(:,1), X_2d(:,2), zeros(N,1),'bo'); 
% grid MINOR
% shg
% 
% 
% 
% % strike =35;
% % dip=40;
% % nhypos =N;
% % 
% % xs_now = xs - mean(xs); % removing the mean for the rotation to be about origins. 
% % ys_now = ys - mean(ys);
% % zs_now = zs - mean(zs);
% % 
% % R=[xs_now ; ys_now ;zs_now];
% % 
% % % rotate into strike direction
% % %Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
% % Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];
% % 
% % Rstrike=Dstrike*R;
% % 
% % rxp(1:nhypos) = Rstrike(1,1:nhypos);
% % ryp(1:nhypos) = Rstrike(2,1:nhypos);
% % rzp(1:nhypos) = Rstrike(3,1:nhypos);
% % 
% % % rotate into dip direction
% % Rstrike(1,1:nhypos) = Rstrike(1,1:nhypos) - mean(Rstrike(1,1:nhypos));
% % Rstrike(2,1:nhypos) = Rstrike(2,1:nhypos) - mean(Rstrike(2,1:nhypos));
% % Rstrike(3,1:nhypos) = Rstrike(3,1:nhypos) - mean(Rstrike(3,1:nhypos));
% % 
% % % Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
% % %Ddip=[ cos(dip) 0 sin(dip); 0 1  0 ; -sin(dip) 0 cos(dip)];
% % Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];
% % 
% % Rdip=Ddip*Rstrike;
% % 
% % Rdip(1,1:nhypos) = Rdip(1,1:nhypos) + mean(Rstrike(1,1:nhypos));
% % Rdip(2,1:nhypos) = Rdip(2,1:nhypos) + mean(Rstrike(2,1:nhypos));
% % Rdip(3,1:nhypos) = Rdip(3,1:nhypos) + mean(Rstrike(3,1:nhypos));
% % 
% % rxp(1:nhypos) = Rdip(1,1:nhypos);
% % ryp(1:nhypos) = Rdip(2,1:nhypos);
% % rzp(1:nhypos) = Rdip(3,1:nhypos);
% % 
% % % Take the depth to have a minimu of zero. Thi does not overright the
% % % original dataset.
% % rzp = rzp-min(rzp);
% % 
% % % Concatenate these result with the original dataet to kep track of each
% % % point. Sort the data based on the new depths.
% % xs_2d = rxp';
% % ys_2d = ryp';
% % zs_2d = rzp';
% % 
% % X = [xs' ys' zs'];
% % 
% %         
% %         
% % figure
% % subplot(1,2,1)
% % plot3(X(:,1), X(:,2), X(:,3),'bo'); %
% % subplot(1,2,2)
% % plot3(xs_2d, ys_2d, zs_2d,'bo'); view(-90,0)
% % grid MINOR
% % shg      
% %         
% %         
%         
%         
% 
% % 
% % 
% % az = -53;el = 31;
% % 
% % V(1,1) = cos(el)*cos(az);
% % V(2,1) = cos(el)*sin(az);
% % V(3,1) = sin(el);
% % 
% % X = [xs' ys' zs'];
% % 
% % X_2d = X - (X*V)*V';
% % 
% % figure
% % subplot(1,2,1)
% % plot3(X(:,1), X(:,2), X(:,3),'bo'); view(az,el)
% % subplot(1,2,2)
% % plot3(X_2d(:,1), X_2d(:,2), zeros(N,1),'bo'); 
% % grid MINOR
% % shg
% 
% 
% 
% 
% % ii = 0;
% % for i=0:5:360
% %     for j=-90:5:90
% %         ii = ii +1;
% %         view(i, j);shg;
% %     end
% % end
% 
% 
