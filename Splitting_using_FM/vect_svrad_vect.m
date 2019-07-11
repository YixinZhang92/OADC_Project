function [svamp_vec, svamp_r] = vect_svrad_vect(m1,m2,m3,m4,m5,m6,pol,vp,vs,dens, r_displ)
%pol = sv_pol;
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
%eps = 1;
con1=(10^5)/(4.0*pi*dens*vs);
hr1=con1.*eps.*p.*eb;                              %vertical strike-slip
hr2=-con1*(eb.*eb - p.*p);                         %vertical dip-slip
hr3=-con1*3.0.*eps.*p.*eb;                         %clvd
% 
% hr1=con1*p.*eb;                                      %vertical strike-slip
% hr2=-con1.*eps.*(eb.*eb - p.*p);                     %vertical dip-slip
% hr3=-con1*3.0.*p.*eb;                                %clvd
%
sa=sin(a);
ca=cos(a);
s2a=sin(2.*a);
c2a=cos(2.*a);
%
%  zero the p_sph array
svrad=zeros(aleng,length(m1));
for i = 1:aleng
      svrad(i,:)= ...
        + (0.5*(m2-m1)*c2a(i) - m4*s2a(i))*hr1(i) ...
        + (m5*ca(i) + m6*sa(i))*hr2(i) ...
        + ((m1+m2-2*m3)/6)*hr3(i);
end
%    
if r_displ == 1

    svamp_z = -p.*svrad; % upward positive
    %svamp_r =  eb.*eps.*svrad; % upward positive %.*eps
    
    %svamp_vec = sign(svamp_r).*sqrt(svamp_z.^2 + svamp_r.^2); %
    
     svamp_r = -eb.*eps.*svrad; % upward positive %.*eps
    
    svamp_vec = sign(svamp_r).*sqrt(svamp_z.^2 + svamp_r.^2); %
    
else
    svamp_vec = svrad;
    svamp_r = svrad;
end
  