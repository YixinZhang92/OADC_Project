function [x,chiscaled] = create_theor_chi_dist(N,df)
%clear all; close all; clc;

%N = 4442;
x = 0:(1/N):4;
%df = 3;


%chisq = sqrt((a*var).^((n/2)-1).* exp(-a/2)/(gamma(n/2)*((2*var)^(n/2))));

chi = (x.^(df-1).*exp(-(x.^2)/2))/(2^((df/2)-1)*gamma(df/2));

nor = 0.9*sqrt(N)/((max(chi))^2)*ones(1,length(chi));
chiscaled = nor.*chi;
 
mean = sqrt(2)*(gamma((df+1)/2)/gamma(df/2));
%trapz(x,chiscaled)
%trapz(chi)

% figure;
% 
% plot(x,chiscaled); shg
% ylim([0 150]); grid MINOR

