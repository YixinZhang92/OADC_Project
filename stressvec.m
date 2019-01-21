function svec=stressvec(m);
%
%   function svec=stressvec(m)
%
%   Calculate the P, T, and I vectors for a moment tensor
%
%   Input:  m = moment tensor array
%           M11=m(1)
%           M22=m(2)
%           M33=m(3)
%           M12=m(4)
%           M13=m(5)
%           M23=m(6)
%
%   Output:  svec = 3x2 array of unit vectors giving the azimuth and
%                   incidence angle of the particular stress axis
%            svec(1,*) = P axis
%            svec(2,*) = T axis
%            svec(3,*) = I axis
%*************************************************************************
%   decompose the moment tensor into eigenvectors and eigenvalues
%   first create the moment tensor 3x3 matrix
c(1,1)=m(1);
c(1,2)=m(4);
c(1,3)=m(5);
c(2,1)=m(4);
c(2,2)=m(2);
c(2,3)=m(6);
c(3,1)=m(5);
c(3,2)=m(6);
c(3,3)=m(3);
%
[v1,d1]=eig(c);						    %	eigenvalue/eigenvectors
%d1
%v1
[v,d]=order(v1,d1);					    %	order the eigenvalues and eigenvectors
%d
%v
ang1=atan2(v(2,1),v(1,1)) * 180/pi;	    % azimuth for first  eigenvalue
ang2=atan2(v(2,2),v(1,2)) * 180/pi;     % azimuth for second eigenvalue
ang3=atan2(v(2,3),v(1,3)) * 180/pi;	    % azimuth for third  eigenvalue
%
vang1=acos(v(3,1))* 180/pi;	        %angle from the vertical for first eigenvalue
vang2=acos(v(3,2))* 180/pi;         %angle from the vertical for second eigenvalue
vang3=acos(v(3,3))* 180/pi;         %angle from the vertical for second eigenvalue
%
%  change to lower hemisphere
%  Note: Matlab uses a right handed cartesian coordinate system such that
%        z is upwards.  The moment tensor cooridinate system is x - north,
%        y - east, and z downward.  Thus, azimuth (clockwise from north is
%        computed correctly in the x/y plane as well as the angle of incidence.
%        An angle of incidence of >= 90 degrees requires that the azimuth
%        be reflected through the origin (add 180 degrees)as well as the
%        complement of the vertical angle.
%
if vang1 > 90.0;
    vang1=180.0-vang1;
    ang1=ang1+180.;
end;
if vang2 > 90.0;
    vang2=180.0-vang2;
    ang2=ang2+180.;
end;
if vang3 > 90.0;
    vang3=180.0-vang3;
    ang3=ang3+180.;
end;
%
%  load svec array
%
svec(1,1)=ang1;
svec(1,2)=vang1;
svec(2,1)=ang2;
svec(2,2)=vang2;
svec(3,1)=ang3;
svec(3,2)=vang3;


return;
