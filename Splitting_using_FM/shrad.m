function [shamp]=shrad(az,inc,vp,vs,dens,eps,m, t_displ)
%
%   function [shamp]=shrad(az,inc,vp,vs,dens,eps,m)
%
%   calculate the SH wave radiation pattern
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
%   Output:   shamp = array of amplitudes
%             shamp(k,j) where k=length(inc); j=length(az);
%
%**************************************************************************
%
%  zero the sh_sph array
%
%sh_sph=zeros(300,1000);
%
con=pi/180.;
p=sin(inc*con)/vs;              %ray parameter array
a=az*con;                       %azimuth array
eb=sqrt(1/(vs*vs) - p.^2);      %vertical slowness array
%
ileng=length(inc);
aleng=length(az);
%
%  compute spherical SV wave source Green's functions
%
%  The resulting amplitude is in cm at 1 km distance
%  from a source with a seismic moment of 10^25 dyn cm.
%
con1=(10^5)/(4.0*pi*dens*vs);
hr1=con1*p/(vs);                              %vertical strike-slip
hr2=-con1.*eps.*eb/(vs);                        %vertical dip-slip
%
%   compute displacements
%
for k=1:ileng
    for j=1:aleng
        sh_sph(k,j)= ...
        + (0.5*(m(1)-m(2))*sin(2*a(j)) - m(4)*cos(2*a(j)))*hr1(k) ...
        + (m(6)*cos(a(j)) - m(5)*sin(a(j)))*hr2(k);
    end
end
shamp=sh_sph;

if t_displ == 1
    shamp = p'.*sh_sph;
end

return;