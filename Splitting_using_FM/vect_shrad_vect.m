function shrad = vect_shrad_vect(m1,m2,m3,m4,m5,m6,pol,vp,vs,dens, t_displ)
%pol = sh_pol;

%  Compute P wave radiation pattern at the azimuth and
%  incidence angle at seismic stations
az = pol(:,2); inc = pol(:,3);
eps = pol(:,6);
%eps = 1;
%
con=pi/180.;
p=sin(inc*con)/vs;              %ray parameter array
a=az*con;                       %azimuth array
eb=sqrt(1/(vs*vs) - p.^2);      %vertical slowness array
%
%ileng=length(inc);
aleng=length(az);
%
%  compute spherical P wave source Green's functions
%
%  The resulting amplitude is in cm at 1 km distance
%  from a source with a seismic moment of 10^25 dyn cm.
%
% Z-axis is positive downward.
con1=(10^5)/(4.0*pi*dens*vs);
hr1=con1*p/(vs);                              %vertical strike-slip
hr2=-con1.*eps.*eb/(vs);                        %vertical dip-slip
%
sa=sin(a);
ca=cos(a);
s2a=sin(2.*a);
c2a=cos(2.*a);
%
%  zero the p_sph array
shrad=zeros(aleng,length(m1));
for i = 1:aleng
      shrad(i,:)= ...
        + (0.5*(m1-m2)*s2a(i) - m4*c2a(i))*hr1(i) ...
        + (m6*ca(i) - m5*sa(i))*hr2(i);
end

if t_displ == 1
    shrad = p.*shrad;
end
