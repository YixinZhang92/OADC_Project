function data = create_syn_rnd_hypo_2D(nhypos,mu,sigma,outfile)
% clear all; close all; clc
% mu= 0;
% sigma=1;
% outfile='rnd_hypos_2D.txt';
% nhypos = 1000;


% Create some random data
s = [mu mu];
x = randn(1000,1);
y1 = normrnd(s(1).*x,sigma);
y2 = normrnd(s(2).*x,sigma);

y1 = y1(randperm(length(y1)));
y2 = y2(randperm(length(y2)));

data = [y1 y2];

%plot(y1,y2,'ro'); shg

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:nhypos
    
    fprintf(fid,'%12.5f %12.5f\n',[y1(kk) y2(kk)]);
    
end