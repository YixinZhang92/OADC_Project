function prad = vect_prad_vect(m1,m2,m3,m4,m5,m6,pol,vp,vs,dens, z_displ)
%pol = p_pol;
%  Compute P wave radiation pattern at the azimuth and
%  incidence angle at seismic stations
az = pol(:,2); inc = pol(:,3);
eps = pol(:,6);
%eps = 1;

con=pi/180.;
p=sin(inc*con)/vp;              %ray parameter array
a=az*con;                       %azimuth array
eba=sqrt(1/(vp*vp) - p.^2);     %vertical slowness array
%
%ileng=length(inc);
aleng=length(az);
%
%  compute spherical P wave source Green's functions
%
%  The resulting amplitude is in cm at 1 km distance
%  from a source with a seismic moment of 10^25 dyn cm.
%
con1=(10^5)/(4.0*pi*dens*vp);
hr0=con1*1.0/(vp*vp);                   %isotropic source
hr1=con1*p.*p;                          %vertical strike-slip
hr2=-2*con1.*eps.*p.*eba;               %vertical dip-slip
hr3=-con1*(p.*p -2.*eba.*eba);          %clvd
%
sa=sin(a);
ca=cos(a);
s2a=sin(2.*a);
c2a=cos(2.*a);
%
%  zero the p_sph array
prad=zeros(aleng,length(m1));

for i = 1:aleng
    prad(i,:) = ((m1+m2+m3)/3)*hr0 ...
        + (0.5*(m2-m1)*c2a(i) - m4*s2a(i))*hr1(i) ...
        + (m5*ca(i) + m6*sa(i))*hr2(i) ...
        + ((m1+m2-2*m3)/6)*hr3(i);
end

if z_displ == 1
    prad = eba.*prad; % positive upward  .*eps
end