% function [c_all,icontour_all, numstart_all, numfin_all] = ...
%     calc_P_nodal_planes(strike,dip,rake,vp,vs,dens, z_displ, r_displ, t_displ)
% % This function determines the P-nodal planes for the strike, dip and rake.
% % The strike, dip and rake can be a vector but of the same size.
% % 
% % Outputs: 
% %            c will contain no of roles = 2 x no of FM
% %            The other outputs will have no roles = no of FM
% %

clear all; close all; clc;


strike = 30; dip =30; rake=30;



N = length(strike);
%
c_all = zeros(2*N,5000);
icontour_all = zeros(N,1);
numstart_all = zeros(N,5);
numfin_all = zeros(N,5);
%
vp=6.0; vs=3.5; dens=2.7; eps=+1.0; z_displ=1;

for no_of_FM = 1:N
        
        %   calculate the moment tensor for a point dislocation
        m=dismom(strike(no_of_FM),dip(no_of_FM),rake(no_of_FM));
        %
        %  Compute P wave radiation pattern
        az=linspace(0,360,1000);
        inc=linspace(0,90,300);
        pamp=prad(az,inc,vp,vs,dens,eps,m, z_displ);
        %
        %  find the zero amplitude contor, if it exists
        c=contourc(az,inc,pamp,[0 0]);
        %
        %  find all of the contours
        [ix,iy]=size(c);
        icontour=0;
        num=0;
        numstart=0;
        numfin=0;
        if iy > 0
        num1=c(2,1);
        if (num1 > 1)
            icontour=1;
            num(1)=num1;
            numstart(1)=2;
            numfin(1)=numstart(1)+num1-1;
        else
            icontour=0;
        end
        size(numfin);
        istop=0;
        ik=0;
        if (icontour > 0)
            test=numfin(1);
            while (test < iy)  
                ik=ik+1;
                icontour=icontour+1;
                num(ik+1)=c(2,numfin(ik)+1);
                numstart(ik+1)=numfin(ik)+2;
                numfin(ik+1)=numstart(ik+1)+num(ik+1)-1;
                test=numfin(ik+1);
            end
        end
        end
 
        [nsr, nsc] = size(numstart);
        [nfr, nfc] = size(numfin);
        icontour_all(no_of_FM,:) = icontour;
        numstart_all(no_of_FM,1:nsc) = numstart;
        numfin_all(no_of_FM,1:nfc) = numfin;
        
        [mmm,nnn]=size(c);
        a = 2*no_of_FM - 1;        b = 2*no_of_FM;
  
        c_all(a:b,1:nnn) = c;       
end

line_color='r';
foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,line_color)












function [m]=dismom(strike,dip,rake)
%
%   function dismom(strike,dip,rake)
%
%   compute the moment tensor for a point dislocation
%   
%
%  Input:   strike = strike of fault plane (clockwise from north, dip to the
%                    right)
%           dip    = fault dip from horizontal
%           rake   = fault rake (slip vector on fault showing slip of upper
%                    plane
%           all angles in degrees
%
%  Output:  M11=m(1)
%           M22=m(2)
%           M33=m(3)
%           M12=m(4)
%           M13=m(5)
%           M23=m(6)
%*************************************************************************
con=pi/180.;
s=strike*con;
d=dip*con;
r=rake*con;
%
m(1)=sin(s)*sin(s)*sin(r)*sin(2*d) + sin(2*s)*cos(r)*sin(d);
m(2)=cos(s)*cos(s)*sin(r)*sin(2*d) - sin(2*s)*cos(r)*sin(d);
m(3)=-sin(r)*sin(2*d);
m(4)=-cos(2*s)*cos(r)*sin(d) - 0.5*sin(2*s)*sin(r)*sin(2*d);
m(5)=cos(s)*cos(r)*cos(d) + sin(s)*sin(r)*cos(2*d);
m(6)=sin(s)*cos(r)*cos(d) - cos(s)*sin(r)*cos(2*d);
%
return;
end


function [pamp]=prad(az,inc,vp,vs,dens,eps,m, z_displ)
%
%   function [pamp]=prad(az,inc,vp,vs,dens,eps,m)
%
%   calculate the radial pwave radiation pattern
%   given an azimuth and incidence angle
%
%   Input:   az = array of azimuths (degrees)
%            inc = array of incidence angles (degrees)
%            eps = array of ray direction values
%                = +1, downgoing ray
%                = -1, upgoing ray
%            vp = source p wave velocity
%            vs = source s wave velocity
%            dens = source density
%            m = moment tensor elements for source (e.g., computed using
%                                                   dismom);
%   Output:   pamp = array of amplitudes
%             pamp(k,j) where k=length(inc); j=length(az);
%
%**************************************************************************
%
%  zero the p_sph array
%
%p_sph=zeros(300,1000);
%
con=pi/180.;
p=sin(inc*con)/vp;              %ray parameter array
a=az*con;                       %azimuth array
eba=real(sqrt(1/(vp*vp) - p.^2));     %vertical slowness array
%save eba.mat eba
%
ileng=length(inc);
aleng=length(az);
%
%  compute spherical P wave source Green's functions
%
%  The resulting amplitude is in cm at 1 km distance
%  from a source with a seismic moment of 10^25 dyn cm.
%
con1=(10^5)/(4.0*pi*dens*vp);
hr0=con1*1.0/(vp*vp);                       %isotropic source
hr1=con1*p.*p;                              %vertical strike-slip
hr2=-2*con1.*eps.*p.*eba;                   %vertical dip-slip
hr3=-con1*(p.*p -2.*eba.*eba);              %clvd
%
%   compute displacements
%
for k=1:ileng
    for j=1:aleng
        p_sph(k,j)=((m(1)+m(2)+m(3))/3)*hr0 ...
        + (0.5*(m(2)-m(1))*cos(2*a(j)) - m(4)*sin(2*a(j)))*hr1(k) ...
        + (m(5)*cos(a(j)) + m(6)*sin(a(j)))*hr2(k) ...
        + ((m(1)+m(2)-2*m(3))/6)*hr3(k);
    end
end
pamp=p_sph;

if z_displ == 1
    pamp = eba'.*p_sph; % positive upward  .*eps
end


return;
end


function foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,line_color)
%
% numstart_all
% numfin_all
%   produce a figure of the zero amplitude contour
%   on an equal area, lower hemisphere, projection
%
%   Input:  c = contour array of zero amplitude contour
%           icontour = number of contours in plot
%           num(icontour)=array of data numbers per contour
%           numstart(icontour)=array of starting positions
%           numfin(icontour)=array of finishing positions for
%                            each contour
%           p_pol = array of P wave polarity data
%                  arranged in polarity, azimuth, incidence
%                  angle 3-tuples
%                  polarity varies as -1,0,+1
%           data_flag = 1 for plot P-wave polarity
%                       0 otherwise
%*********************************************************
%
%   create the outside circle of the projection
%
con=pi/180.;
radius=1.0;
a=linspace(0,360,1000)*con;
yc=radius*cos(a);
xc=radius*sin(a);


% This function updates the azimuth at each station based on the direction
% of wave propagation. 
% az = az + 180 for upgoing waves
% az = az       for downgoing waves
%p_pol = update_az(p_pol);

% figure;
% figure('color','w','name','P Wave Nodal Surfaces');
% hold on;
% daspect([1 1 1]);
% axis off;
% plot(xc,yc,'linewidth',1.5,'color',[0.5,0.5,0.5]);

%figure;
figure('color','w','name','P Wave Nodal Surfaces');
%hold on;
daspect([1 1 1]);
%axis off;
plot(xc,yc,'linewidth',1.5,'color',[0.5,0.5,0.5]);


[no_of_FM,~] = size(icontour_all);

for j = 1:no_of_FM
    icontour = icontour_all(j,:);
    a = 2*j - 1; b = 2*j;
        %j
    for m=1:icontour
        %m
        c = c_all(a:b,:); c( :, ~any(c,1) ) = []; % Removing columns containing zeros
        numstart = numstart_all(j,:); numstart( :, ~any(numstart,1) ) = [];
        numfin = numfin_all(j,:); numfin( :, ~any(numfin,1) ) = [];
        
%         if j == 1
%         save c.mat c;
%         save numstart.mat numstart;
%         save numfin.mat numfin;
%         end
        
        azn=c(1,numstart(m):numfin(m));
        incn=c(2,numstart(m):numfin(m));
        r=sqrt(2)*radius*sin(incn*con/2);
        ycnod=r.*cos(azn*con);
        xcnod=r.*sin(azn*con);
        plot(xcnod,ycnod,'linewidth',1.0,'color', line_color);        
    end
end
%
%   plot N,E,S,W and center tick marks
%
for it=1:4
    az=90*(it-1)*con;
    r1=radius*0.97;
    r2=radius*1.03;
    xt=radius*1.1*sin(az);
    yt=radius*1.1*cos(az);
    plot([r1*sin(az) r2*sin(az)],[r1*cos(az) r2*cos(az)]);
    if it == 1
        text(xt,yt,'N','horizontalalignment','center', 'FontSize',14);
    end
    if it == 2
        text(xt,yt,'E','horizontalalignment','center', 'FontSize',14);
    end
    if it == 3
        text(xt,yt,'S','horizontalalignment','center', 'FontSize',14);
    end
    if it == 4
        text(xt,yt,'W','horizontalalignment','center', 'FontSize',14);
    end   
end
plot([-0.03 0.03],[0 0]);
plot([0 0],[-0.03 0.03]);
%
title('P Wave Nodal Surfaces', 'FontSize', 20);

hold off;
return;
end





