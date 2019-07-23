function shfocsphere(nmech,mom,stype,vp,vs,dens,sh_pol,nodeflag,axflag,dataflag, t_displ)
%
%   function shfocsphere(m,stype,vp,vs,dens,sh_pol,nodeflag,axflag,dataflag)
%
%   Plot the lower hemisphere equal area projection of the SV wave
%   focal sphere.
%
%   Input:  nmech = number of focal mechanisms to be plotted.
%
%           mom = (nmechx6) array of moment tensor elements if stype=0
%           M11=mom(1)
%           M22=mom(2)
%           M33=mom(3)
%           M12=mom(4)
%           M13=mom(5)
%           M23=mom(6)
%                          -OR-
%           mom= (nmechx3) array containing dislocation parameters if stype=1
%           mom(1)= fault strike
%           mom(2)= fault dip
%           mom(3)= fault rake
%                using the Langston and Helmberger(1975) conventions
%                Multiple nodal planes can be plotted (nmech>1)
%
%           vp, vs, dens = source P wave, S wave velocities (km/s)
%                          and density (gm/cc)
%
%           sh_pol = (nx3) array
%                  arranged in polarity, azimuth, incidence angle 3-tuples
%                  polarity varies as -1,0,+1 (0 is nodal)
%                  
%
%           nodeflag = 1, plot the nodal surfaces
%                    = 0, do not plot the nodal surfaces
%
%           axflag = 1, plot the positions of the P, T, and I axes
%                  = 0, do not plot the axes
%
%           dataflag = 1, plot the positions of the data polarities
%                    = 0, do not plot the data
%
%   Ouput:  A figure displaying mechanisms nodal planes, stress axes, and
%           data.
%
%   C.A. Langston 9/11/03
%**************************************************************************
%
%    Prepare the plot figure
%
%**************************************************************************
%   create the outside circle of the projection
%
con=pi/180.;
radius=1.0;
a=linspace(0,360,1000)*con;
yc=radius*cos(a);
xc=radius*sin(a);
%
% This function updates the azimuth at each station based on the direction
% of wave propagation. 
% az = az + 180 for upgoing waves
% az = az       for downgoing waves
%p_pol = update_az(p_pol);

I = (sh_pol(:,6) == -1);
sh_pol(I,2) = sh_pol(I,2)+180;

% converting P-wave amplitudes to polarities
sh_pol(:,1) = sh_pol(:,1)./abs(sh_pol(:,1));

%hshfoc=figure('color','w','name','SH Wave Nodal Surfaces');
hold on;
daspect([1 1 1]);
%axis off;
plot(xc,yc,'linewidth',1.5,'color',[0.5,0.5,0.5]);
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
    if it == 1;
        text(xt,yt,'N','horizontalalignment','center');
    end;
    if it == 2;
        text(xt,yt,'E','horizontalalignment','center');
    end;
    if it == 3;
        text(xt,yt,'S','horizontalalignment','center');
    end;
    if it == 4;
        text(xt,yt,'W','horizontalalignment','center');
    end;   
end;
plot([-0.03 0.03],[0 0]);
plot([0 0],[-0.03 0.03]);
%
title('SH Wave Nodal Surfaces');
%
%*************************************************************************
%   Plot the data, if wanted
%*************************************************************************
%
if dataflag == 1;
%  Put polarity info onto plot
plot(0.9,1.2,'marker','o','markeredgecolor','k','markerfacecolor','k');
text(1.0,1.2,'clockwise');
plot(0.9,1.05,'marker','x','markeredgecolor','k','markerfacecolor','k');
text(1.0,1.05,'Nodal');
plot(0.9,0.90,'marker','o','markeredgecolor','k','markerfacecolor','w');
text(1.0,0.90,'counterclockwise');
%
%   plot wave polarities
%
[pix,piy]=size(sh_pol);
ik=1;
while ik <= pix;
    pol=sh_pol(ik,1);
    az=sh_pol(ik,2);
    inc=sh_pol(ik,3);
    r=sqrt(2)*radius*sin(inc*con/2);
    xp=r*sin(az*con);
    yp=r*cos(az*con);
    if pol == -1.0;
        mark='o';
        fcol='w';
    end;
    if pol == 0.0;
        mark='x';
        fcol='k';
    end;
    if pol == 1.0;
        mark='o';
        fcol='k';
    end;
    plot(xp,yp,'marker',mark,'markeredgecolor','k','markerfacecolor',fcol);
    ik=ik+1;
end;
end;
%
%*************************************************************************
%   Now for the mechanism loop (will assume either/both axes and nodal 
%                               surfaces are wanted)
%*************************************************************************
%
for imech=1:nmech;
    if stype == 1;
        m=dismom(mom(imech,1),mom(imech,2),mom(imech,3));    % dislocation model
    else;
        m(1:6)=mom(imech,1:6);
    end;
    %
    %   Decide on and plot stress axes for this mechanism
    %
    if axflag == 1;
        %
        %   calculate the P, T, and I vectors for this moment tensor
        %
        svec=stressvec(m);
        %
        %   Plot P stress axes
        %
        rp=sqrt(2.0)*radius*sin(svec(1,2)*con/2.);
        yct=rp*cos(svec(1,1)*con);
        xct=rp*sin(svec(1,1)*con);
        text(xct,yct,'P','horizontalalignment','center');
        %
        %   Plot I stress axes
        %
        rp=sqrt(2.0)*radius*sin(svec(2,2)*con/2.);
        yct=rp*cos(svec(2,1)*con);
        xct=rp*sin(svec(2,1)*con);
        text(xct,yct,'I','horizontalalignment','center');
        %
        %   Plot T stress axes
        %
        rp=sqrt(2.0)*radius*sin(svec(3,2)*con/2.);
        yct=rp*cos(svec(3,1)*con);
        xct=rp*sin(svec(3,1)*con);
        text(xct,yct,'T','horizontalalignment','center');
    end;
    %
    %   Decide on and plot nodal planes for this mechanism
    %
    if nodeflag == 1;
        %
        %  Compute P wave radiation pattern
        %
        eps=+1.0;
        az=linspace(0,360,1000);
        inc=linspace(0,90,300);
        shamp=shrad(az,inc,vp,vs,dens,eps,m, t_displ);
        %
        %  find the zero amplitude contor, if it exists
        %
        c=contourc(az,inc,shamp,[0 0]);
        %  find all of the contours
        [ix,iy]=size(c);
        icontour=0;
        num=0;
        numstart=0;
        numfin=0;
        if iy > 0;
            num1=c(2,1);
            if (num1 > 1);
                icontour=1;
                num(1)=num1;
                numstart(1)=2;
                numfin(1)=numstart(1)+num1-1;
            else;
                icontour=0;
            end;
            size(numfin);
            istop=0;
            ik=0;
        if (icontour > 0);
            test=numfin(1);
            while (test < iy);  
                ik=ik+1;
                icontour=icontour+1;
                num(ik+1)=c(2,numfin(ik)+1);
                numstart(ik+1)=numfin(ik)+2;
                numfin(ik+1)=numstart(ik+1)+num(ik+1)-1;
                test=numfin(ik+1);
             end;
        end;
        end;
        %
        %    load and plot each nodal contour
        %
        for m=1:icontour;
            azn=c(1,numstart(m):numfin(m));
            incn=c(2,numstart(m):numfin(m));
            r=sqrt(2)*radius*sin(incn*con/2);
            ycnod=r.*cos(azn*con);
            xcnod=r.*sin(azn*con);
            plot(xcnod,ycnod,'linewidth',1.0,'color','r');
        end;
    end;
end;
%
%   All done with plot
%
%hold off;
return;