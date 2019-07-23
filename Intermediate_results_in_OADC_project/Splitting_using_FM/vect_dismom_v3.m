function [m1,m2,m3,m4,m5,m6]=vect_dismom_v3(strike,dip,rake)
%
%   function dismom_v3(strike,dip,rake)
%
%   compute the moment tensor for a point dislocation - vectorize
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
%
% Seun changed the sign in m2 to negative from positive and flattens the moment tensor components.
%*************************************************************************
con=pi/180.;
s=strike*con;
d=dip*con;
r=rake*con;

ns=length(strike);
nr=length(rake);
nd=length(dip);

m1(1:nd,1:ns,1:nr)=0.0;
m2(1:nd,1:ns,1:nr)=0.0;
m3(1:nd,1:ns,1:nr)=0.0;
m4(1:nd,1:ns,1:nr)=0.0;
m5(1:nd,1:ns,1:nr)=0.0;
m6(1:nd,1:ns,1:nr)=0.0;

ss=sin(s);
cs=cos(s);
s2s=sin(2.*s);
c2s=cos(2.*s);
sr=sin(r);
cr=cos(r);
sd=sin(d);
cd=cos(d);
s2d=sin(2.*d);
c2d=cos(2.*d);

A1=(ss.*ss)'*sr;
B1=s2s'*cr;

A2=(cs.*cs)'*sr;

s1(1:ns)=1.0;
A3=-s1'*sr;

A4=-c2s'*cr;
B4=-0.5.*s2s'*sr;

A5=cs'*cr;
B5=ss'*sr;

A6=ss'*cr;
B6=-cs'*sr;

for k=1:nd
    m1(k,1:ns,1:nr)=s2d(k).*A1(1:ns,1:nr) + sd(k).*B1(1:ns,1:nr);
    m2(k,1:ns,1:nr)=s2d(k).*A2(1:ns,1:nr) - sd(k).*B1(1:ns,1:nr);
    m3(k,1:ns,1:nr)=s2d(k).*A3(1:ns,1:nr);
    m4(k,1:ns,1:nr)=sd(k).*A4(1:ns,1:nr) + s2d(k).*B4(1:ns,1:nr);
    m5(k,1:ns,1:nr)=cd(k).*A5(1:ns,1:nr) + c2d(k).*B5(1:ns,1:nr);
    m6(k,1:ns,1:nr)=cd(k).*A6(1:ns,1:nr) + c2d(k).*B6(1:ns,1:nr);
end
    

% I flattens the moment tensor components. It is easier to determine the 
% displacements 
m1 = vect_func_flattens(m1);
m2 = vect_func_flattens(m2);
m3 = vect_func_flattens(m3);
m4 = vect_func_flattens(m4);
m5 = vect_func_flattens(m5);
m6 = vect_func_flattens(m6);
