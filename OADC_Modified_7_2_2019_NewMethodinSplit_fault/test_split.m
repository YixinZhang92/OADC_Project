figure
plot(xv(1,1),yv(1,1),'ro'); hold on;
plot(xv(1,2),yv(1,2),'bo'); hold on;
plot(xv(1,3),yv(1,3),'go'); hold on;
plot(xv(1,4),yv(1,4),'ko'); hold on;

plot(xv1,yv1,'ro'); hold on;
plot(xv2,yv2,'bo'); hold on;
plot(xv3,yv3,'go'); hold on;
plot(xv4,yv4,'ko'); hold on;


xm1 = (xv(1,1) + xv(1,2))/2;
xm2 = (xv(1,3) + xv(1,4))/2;

ym1 = (yv(1,1) + yv(1,2))/2;
ym2 = (yv(1,3) + yv(1,4))/2;

plot(xm1,ym1,'ro'); hold on;
plot(xm2,ym2,'ro'); hold on;


xr1 = (3*xv(1,1) + xv(1,2) + xv(1,3) +3*xv(1,4))/8;
xr2 = (xv(1,1) + 3*xv(1,2) + 3*xv(1,3) +xv(1,4))/8;


yr1 = (3*yv(1,1) + yv(1,2) + yv(1,3) +3*yv(1,4))/8;
yr2 = (yv(1,1) + 3*yv(1,2) + 3*yv(1,3) +yv(1,4))/8;

plot(xr1,yr1,'go'); hold on;
plot(xr2,yr2,'ko'); hold on;