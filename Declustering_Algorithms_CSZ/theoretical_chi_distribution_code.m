clear all; close all; clc;

N = 4442;
x = 0:(1/N):4;
k = 3;


%chisq = sqrt((a*var).^((n/2)-1).* exp(-a/2)/(gamma(n/2)*((2*var)^(n/2))));

chi = (x.^(k-1).*exp(-(x.^2)/2))/(2^((k/2)-1)*gamma(k/2));

nor = 0.9*sqrt(N)/((max(chi))^2)*ones(1,length(chi));
chiscaled = nor.*chi;
 
mean = sqrt(2)*(gamma((k+1)/2)/gamma(k/2));
%trapz(x,chiscaled)
%trapz(chi)

figure;

plot(x,chiscaled); shg
ylim([0 150]); grid MINOR

