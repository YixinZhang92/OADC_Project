function [svamp_vec, svamp_r]=svrad(az,inc,vp,vs,dens,eps,m, r_displ)
%
%   function [svamp]=prad(az,inc,vp,vs,dens,eps,m)
%
%   calculate the latitudinal SV wave radiation pattern
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
%   Output:   svamp = array of amplitudes
%             svamp(k,j) where k=length(inc); j=length(az);
%
%**************************************************************************
%
%  zero the sv_sph array
%
%sv_sph=zeros(300,1000);
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
%eps = 1;
con1=(10^5)/(4.0*pi*dens*vs);
hr1=con1*eps.*p.*eb;                              %vertical strike-slip
hr2=-con1*(eb.*eb - p.*p);                        %vertical dip-slip
hr3=-con1*3.0*eps*p.*eb;                          %clvd
%
% hr1=con1*p.*eb;                                     %vertical strike-slip
% hr2=-con1*eps.*(eb.*eb - p.*p);                     %vertical dip-slip
% hr3=-con1*3.0*p.*eb;                                %clvd
%
%   compute displacements
%
for k=1:ileng
    for j=1:aleng
        sv_sph(k,j)= ...
        + (0.5*(m(2)-m(1))*cos(2*a(j)) - m(4)*sin(2*a(j)))*hr1(k) ...
        + (m(5)*cos(a(j)) + m(6)*sin(a(j)))*hr2(k) ...
        + ((m(1)+m(2)-2*m(3))/6)*hr3(k);
    end
end
%svamp=sv_sph;

if r_displ == 1
    
    svamp_z = -p'.*sv_sph; % upward positive
    %svamp_r =  eb'.*eps.*sv_sph; % upward positive %
    
    %svamp_vec = sign(svamp_r).* sqrt(svamp_z.^2 + svamp_r.^2); %sign(svamp_r)
    
    svamp_r = -eb'.*eps.*sv_sph; % upward positive %
    
    svamp_vec = sign(svamp_r).* sqrt(svamp_z.^2 + svamp_r.^2); %
    
else
    svamp_vec = sv_sph;
    svamp_r = sv_sph;
end

return;