function [alph,beta,gama]=trisol(a,b,c,flag)

% Solve for the angles of a triangle given the length of each side

s=0.5.*(a+b+c);
r=sqrt((s-a).*(s-b).*(s-c)./s);

alph=2.0.*atan2(r,(s-a));
beta=2.0.*atan2(r,(s-b));
gama=2.0.*atan2(r,(s-c));

if flag == 'd'
    con=180./pi;
    alph=alph.*con;
    beta=beta.*con;
    gama=gama.*con;
end;
return;

