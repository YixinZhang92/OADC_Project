clear all; close all; clc;
N = 1000;
% Create theoretical chi distribution
[x,chiscaled] = create_theor_chi_dist(N,1);

% stdn = sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2);
% chiscaled = chiscaled*sqrt(N)/stdn;
chiscaled = chiscaled*100;

figure;
plot(x, chiscaled,'r'); hold on;


% synthetic hypocenters
mu= 0; sigma=1; nhypos = 1000;
outfile='rnd_hypos_3D.txt';
infile = outfile;
[data,R] = create_syn_rnd_hypo_3D(nhypos,mu,sigma,outfile);

RR =  sqrt(data(:,1).^2 + data(:,2).^2 + data(:,3).^2);




[a,b] = hist(RR,50);

stairs(b,a,'b','LineWidth',1.5)


figure
pd1 = createFit_cdf(RR)



df=1;
chi = (RR.^(df-1).*exp(-(RR.^2)/2))/(2^((df/2)-1)*gamma(df/2))*100;

plot(RR,chi,'bo'); hold on;






% [a,b] = hist(chi,50); hold on;
% stairs(b,a,'b','LineWidth',1.5)