function [aux_strike, aux_dip, aux_rake]= auxiliary_fault_plane_based_on_Kumar_and_Zhao(strike, dip, rake)
% ! Program to calculate the strike, dip and rake of the auxiliary fault plane solutions
% ! given the strike, dip and rake of the fault plane.
% Originally by Utpal Kumar, Li Zhao in Fortran format
% Seun Fadugba wrote a MATLAB version of the code 

%clear all; close all; clc;

N=length(strike);
[m,n] = size(strike);
aux_strike=zeros(N,1); aux_dip=zeros(N,1); aux_rake=zeros(N,1);

% strike= 278;
% dip= 68;
% rake= -68;

% The strike of the auxiliary plane is:   50.8362 degrees
% The dip of the auxiliary plane is:   30.7205 degrees
% The rake of the auxiliary plane is: -132.8362 degrees


for i=1:N
    ph1=strike(i);
    del1=dip(i);
    lb1=rake(i);
    %Converting the strike, dip and rake in radians
    ph1=ph1*pi/180;  %!strike
    del1=del1*pi/180;  %      !dip
    lb1=lb1*pi/180   ;  %     !rake

    %!Adaptation to include the negative values of rake
    lbsgn=1;
    if (lb1 < 0)
    lbsgn=-1;
    end


    %! Calculation of Strike, sip and rake of the auxiliary fault plane
    del2=acos(sin(lb1)*sin(del1));  % !dip of auxialiary fault plane
    sinlb2=cos(del1)/sin(del2);
    coslb2=-(sin(del1)*cos(lb1)/sin(del2));
    lb2=acos(coslb2)    ;    %!rake of auxilairy fault plane
    sinph1_ph2=cos(lb1)/sin(del2);
    cosph1_ph2=-1/(tan(del1)*tan(del2));
    ph1_ph2=acos(cosph1_ph2);

    %!Checking for the quadrant of the strike angle
    if (sinph1_ph2 >= 0 && cosph1_ph2 >=0)
        ph1_ph2=ph1_ph2;
    elseif (sinph1_ph2 > 0 && cosph1_ph2 < 0)
        ph1_ph2=ph1_ph2;
    elseif (sinph1_ph2 < 0 && cosph1_ph2 < 0) 
        ph1_ph2=-ph1_ph2;
    elseif (sinph1_ph2 < 0 && cosph1_ph2 > 0) 
        ph1_ph2=-ph1_ph2;
    end

    ph2=ph1-ph1_ph2;        % !strike of auxialiary fault plane

    %! For dip > 90 degrees and less than 180 degrees
    if (del2 > pi/2 && del2 < pi)
    ph2= pi + ph2;
    del2= pi - del2;
    lb2= 2*pi - lb2;
    end

    if (lbsgn < 0)
    lb2 = -(2*pi - lb2);
    end

    %!Adaptation to give the strike value in the range of 0 to 360 degrees
    if (ph2 > 2*pi)
    ph2 = ph2 - 2*pi;
    end


    aux_strike(i) = real(ph2*180/pi);
    aux_dip(i)= real(del2*180/pi);
    aux_rake(i) = real(lb2*180/pi);

    if aux_strike(i) < 0
        aux_strike(i) = aux_strike(i)+360;
    end
end

aux_strike = reshape(aux_strike, [m,n]);
aux_dip = reshape(aux_dip, [m,n]);
aux_rake = reshape(aux_rake, [m,n]);
end


%% Original fortran code

% program auxiliary_fault_plane
% ! Program to calculate the strike, dip and rake of the auxiliary fault plane solutions
% ! given the strike, dip and rake of the fault plane.
% ! Authors: Utpal Kumar, Li Zhao
% implicit none
% real (kind=8):: ph1,del1,lb1,ph2,del2,lb2,sinlb2,coslb2,sinph1_ph2,cosph1_ph2,ph1_ph2,phI,delI,lbI
% real, parameter :: pi = 2*asin(1.0)
% integer :: n,lbsgn
% 
% print *, 'Enter strike, dip and rake respectively of the fault plane (in degrees) (e.g., 40,50,60):'
% read *, ph1,del1,lb1
% !ph1=40
% !del1=80
% !lb1=200
% phI=ph1
% delI=del1
% lbI=lb1
% !Converting the strike, dip and rake in radians
% ph1=ph1*pi/180  !strike
% del1=del1*pi/180        !dip
% lb1=lb1*pi/180          !rake
% 
% !Adaptation to include the negative values of rake
% lbsgn=1
% if (lb1 < 0) then
% lbsgn=-1
% end if
% 
% 
% ! Calculation of Strike, sip and rake of the auxiliary fault plane
% del2=acos(sin(lb1)*sin(del1))   !dip of auxialiary fault plane
% sinlb2=cos(del1)/sin(del2)
% coslb2=-(sin(del1)*cos(lb1)/sin(del2))
% lb2=acos(coslb2)        !rake of auxilairy fault plane
% sinph1_ph2=cos(lb1)/sin(del2)
% cosph1_ph2=-1/(tan(del1)*tan(del2))
% ph1_ph2=acos(cosph1_ph2)
% 
% !Checking for the quadrant of the strike angle
% if (sinph1_ph2 >= 0 .and. cosph1_ph2 >=0) then
% ph1_ph2=ph1_ph2
% else if (sinph1_ph2 > 0 .and. cosph1_ph2 < 0) then
% ph1_ph2=ph1_ph2
% else if (sinph1_ph2 < 0 .and. cosph1_ph2 < 0) then
% ph1_ph2=-ph1_ph2
% else if (sinph1_ph2 < 0 .and. cosph1_ph2 > 0) then
% ph1_ph2=-ph1_ph2
% end if
% 
% ph2=ph1-ph1_ph2         !strike of auxialiary fault plane
% 
% ! For dip > 90 degrees and less than 180 degrees
% if (del2 > pi/2 .and. del2 < pi) then
% ph2= pi + ph2
% del2= pi - del2
% lb2= 2*pi - lb2
% end if
% 
% if (lbsgn < 0) then
% lb2 = -(2*pi - lb2)
% end if
% 
% !Adaptation to give the strike value in the range of 0 to 360 degrees
% if (ph2 > 2*pi) then
% ph2 = ph2 - 2*pi
% end if
% 
% 101 format("The strike of the auxiliary plane is: ",f9.4, " degrees")
% 102 format("The dip of the auxiliary plane is: ",f9.4, " degrees")
% 103 format("The rake of the auxiliary plane is: ",f9.4, " degrees")
% print 101, ph2*180/pi
% print 102,  del2*180/pi
% print 103, lb2*180/pi
% open(unit=10,file='plt.dat')
% 10 format("25  25   0 ",  f8.2, f8.2, f8.2, f8.2, f8.2, f8.2," MainFault: ",f7.2,"/",f7.2,"/",f7.2, &    !continuation of line
% " AuxFault: ",f7.2,"/",f7.2,"/",f7.2)
% write(10,10) phI,delI,lbI,ph2*180/pi, del2*180/pi, lb2*180/pi, phI,delI,lbI,ph2*180/pi, del2*180/pi, lb2*180/pi !writing in the file
% end program auxiliary_fault_plane
