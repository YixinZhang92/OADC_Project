X=[xs' ys' zs'];
close all;
% compute the covariance matrix for this cluster
Cxy=cov(X,0);
[V,D]=eig(Cxy);
basis1=V(1:3,2)
basis2=V(1:3,3)
vec_plane=V(1:3,1);

k=1;
% compute the covariance matrix for this cluster
Cxy=cov(X,0);

% compute the eigenvalues and eigenvectors for this cluster
[V,D]=eig(Cxy);
% 
% % Seun changed the way we determine the fault vertex, length and width
% % of the fault.
% % calculate fault plane parameters from the eigen results
% % and calculate the vertices of the fault plane
% 
% % compute vertice locations
% xv(k,1)=min(xs);
% yv(k,1)=max(ys);
% 
% xv(k,2)=min(xs);
% yv(k,2)=min(ys);
% 
% xv(k,3)=max(xs);
% yv(k,3)=min(ys);
% 
% xv(k,4)=max(xs);
% yv(k,4)=max(ys);
% 
% zv(k,:) = (1/V(3,1)) .* ((xb(k)*V(1,1) + yb(k)*V(2,1) + (zb(k)*V(3,1))) ...
%     - (V(1,1)*[xv(k,1) xv(k,2) xv(k,3) xv(k,4)] + V(2,1)*[yv(k,1) yv(k,2) yv(k,3) yv(k,4)] ));
% 
% L(k) = sqrt((xv(k,2)-xv(k,1))^2 + (yv(k,2)-yv(k,1))^2 + (zv(k,2)-zv(k,1))^2);
% W(k) = sqrt((xv(k,4)-xv(k,1))^2 + (yv(k,4)-yv(k,1))^2 + (zv(k,4)-zv(k,1))^2);
% 
% % compute dip and strike
% con=180.0./pi;
% % Dip=con.*acos(V(3,1));
% % Strike=con.*atan2(V(1,1),V(2,1));
% 
% Dip(k)=con.*acos(V(3,1));
% if (V(1,1) >= 0) && (V(2,1) >= 0)
%     Strike(k)=con.*atan(abs(V(1,1))/abs(V(2,1)))-90;
% elseif (V(1,1) >= 0) && (V(2,1) < 0)
%     Strike(k)=180-con.*atan(abs(V(1,1))/abs(V(2,1)))-90;
% elseif (V(1,1) < 0) && (V(2,1) < 0)
%     Strike(k)=180+con.*atan(abs(V(1,1))/abs(V(2,1)))-90;    
% elseif (V(1,1) < 0) && (V(2,1) >= 0)
%     Strike(k)=360-con.*atan(abs(V(1,1))/abs(V(2,1)))-90;    
% end    
% 
% if Strike(k) < 0
%     Strike(k) = Strike(k)+180;
% end
% 
% % save the plane unit normal vector and eigenvalue
% vec_plane(k,1:3)=V(1:3,1);
% %lambda3(k)=sqrt(12.*D(1,1));
% lambda3(k)=sqrt(D(1,1));
% 
% picname='Test L and W 1'; simul_tag = '1';
% datplot(xs,ys,zs,1,xv,yv,zv,picname,simul_tag);
% 

xb = mean(xs);
yb = mean(ys);
zb = mean(zs);

[L(k),W(k),Strike(k),Dip(k),xv(k,:),yv(k,:),zv(k,:)] = fltplane(V,D,xb,yb,zb);

picname='Test L and W 2'; simul_tag = '1';
datplot(xs,ys,zs,1,xv,yv,zv,picname,simul_tag);

% save the plane unit normal vector and eigenvalue
vec_plane(k,1:3)=V(1:3,1);
%lambda3(k)=sqrt(12.*D(1,1));
lambda3(k)=sqrt(D(1,1));






nhypos = length(xs);

for i = 1:nhypos
AP = [0 0 0] + X(i,:);
AB = [0 0 0] - [V(1,3) V(2,3) V(3,3)];
NN(i,:) = [0 0 0] + dot(AP,AB) / dot(AB,AB) * AB;

end

ex_dist = dot([xb yb zb],AB) / dot(AB,AB) * AB;

L = sqrt((min(NN(:,1))-max(NN(:,1)))^2 + (min(NN(:,2))-max(NN(:,2)))^2 +(min(NN(:,3))-max(NN(:,3)))^2);

hold on;
%V = 20*V;
plot3(NN(:,1), NN(:,2), NN(:,3), 'ro')



for i = 1:nhypos
AP = [0 0 0] + X(i,:);
AB = [0 0 0] - [V(1,2) V(2,2) V(3,2)];
NN(i,:) = [0 0 0] + dot(AP,AB) / dot(AB,AB) * AB;

end

W = sqrt((min(NN(:,1))-max(NN(:,1)))^2 + (min(NN(:,2))-max(NN(:,2)))^2 +(min(NN(:,3))-max(NN(:,3)))^2);



hold on;
%V = 20*V;
plot3(NN(:,1), NN(:,2), NN(:,3), 'ro')












hold on;
%plot3(NN(k,1),NN(k,2),NN(k,3),'ko');

hold on;
%V = 20*V;
plot3(NN(:,1), NN(:,2), NN(:,3), 'ro')

hold on;
plot3([0 V(1,3)], [0 V(2,3)], [0 V(3,3)],'g') 


hold on;
plot3([0 V(1,1)], [0 V(2,1)], [0 V(3,1)],'k') 


hold on;
plot3([0 V(1,2)], [0 V(2,2)], [0 V(3,2)],'c') 







R=[NN(:,1) NN(:,2) NN(:,3)];
[n,~] = size(R);


VV = [V(:,3) V(:,2) V(:,1)];


Rn = VV\R';

xvt(1:n)=Rn(1,1:n);
yvt(1:n)=Rn(2,1:n);
zvt(1:n)=Rn(3,1:n);

figure;
plot3(xvt,yvt,zvt,'ro')









V = 20*V;
quiver(xb,yb,V(1,3),V(2,3)); hold on;

quiver(xb,yb,V(1,2),V(2,2)); hold on;






con=pi/180.;

rake = 0;
strike=-Strike(k).*con;
dip=Dip(k).*con;
rake=rake.*con;

xvt = xs-mean(xs);
yvt = ys-mean(ys);
zvt = zs-mean(zs);

R=[xvt; yvt; zvt];
[~,nhypos] = size(R);


        % rotate into strike direction
        %Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
        Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];

        Rstrike=Dstrike*R;

        rxp(1:nhypos) = Rstrike(1,1:nhypos);
        ryp(1:nhypos) = Rstrike(2,1:nhypos);
        rzp(1:nhypos) = Rstrike(3,1:nhypos);

        % rotate into dip direction
        Rstrike(1,1:nhypos) = Rstrike(1,1:nhypos) - mean(Rstrike(1,1:nhypos));
        Rstrike(2,1:nhypos) = Rstrike(2,1:nhypos) - mean(Rstrike(2,1:nhypos));
        Rstrike(3,1:nhypos) = Rstrike(3,1:nhypos) - mean(Rstrike(3,1:nhypos));

         %Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
        %Ddip=[ cos(dip) 0 sin(dip); 0 1  0 ; -sin(dip) 0 cos(dip)];
        Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];

        Rdip=Ddip*Rstrike;

        Rdip(1,1:nhypos) = Rdip(1,1:nhypos) + mean(Rstrike(1,1:nhypos));
        Rdip(2,1:nhypos) = Rdip(2,1:nhypos) + mean(Rstrike(2,1:nhypos));
        Rdip(3,1:nhypos) = Rdip(3,1:nhypos) + mean(Rstrike(3,1:nhypos));

        rxp(1:nhypos) = Rdip(1,1:nhypos);
        ryp(1:nhypos) = Rdip(2,1:nhypos);
        rzp(1:nhypos) = Rdip(3,1:nhypos);




picname='Test L and W 2'; simul_tag = '1';
datplot(rxp,ryp,rzp,1,xv,yv,zv,picname,simul_tag);



































% xvt = xs-mean(xs);
% yvt = ys-mean(ys);
% zvt = zs-mean(zs);

R=[xvt; yvt; zvt];
[~,n] = size(R);


VV = [V(:,3) V(:,2) V(:,1)];


Rn = VV\R;

xvt(1:n)=Rn(1,1:n);
yvt(1:n)=Rn(2,1:n);
zvt(1:n)=Rn(3,1:n);

figure;
plot3(xvt,yvt,zvt,'ro')


%%
con=pi/180.;

rake = 0;
strike=Strike(k).*con;
dip=Dip(k).*con;
rake=rake.*con;

% Fault corner coordinates have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike

% rotate into rake direction
xvt = xs;
yvt = ys;
zvt = zs;

R=[xvt; yvt; zvt];
[~,n] = size(R);

% rotate into strike direction

%Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
%Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];
Dstrike=[ -sin(strike)  -cos(strike) 0 ; cos(strike) -sin(strike) 0 ; 0 0 1];

Rstrike=Dstrike*R;

xvt(1:n)=Rstrike(1,1:n);
yvt(1:n)=Rstrike(2,1:n);
zvt(1:n)=Rstrike(3,1:n);

figure;
plot3(xvt,yvt,zvt,'ro')






% rotate into dip direction

%Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];

Rdip=Ddip*Rstrike;

xvt(1:n)=Rdip(1,1:n);
yvt(1:n)=Rdip(2,1:n);
zvt(1:n)=Rdip(3,1:n);

% picnamet = 'rotate.into.dip';
% 
% datplot(xst,yst,zst,kmint,xvt,yvt,zvt,picnamet,simul_tagt);





figure;
plot3(xvt,yvt,zvt,'ro')