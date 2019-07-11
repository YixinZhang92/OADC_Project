clear all; close all; clc;

a=0:20; % a is scaled by the variance
var = 10;

for n=3%1:10

    %chisq = ((a*var).^((n/2)-1).* exp(-a/2)/(gamma(n/2)*((2*var)^(n/2))));

    chisq = sqrt((a*var).^((n/2)-1).* exp(-a/2)/(gamma(n/2)*((2*var)^(n/2))));
    
    plot(a,chisq); hold on;
    
end
   grid MINOR; shg

