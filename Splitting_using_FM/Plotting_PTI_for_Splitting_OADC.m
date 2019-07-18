clear all; close all; clc;

infile = 'Combined_Dataset_MT_PL.csv';
quality_threshold = 4;

%%
fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[11,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs=data(1:N,1);
ys=data(1:N,2);
zs=-data(1:N,3);
strikes = data(1:N,7);
dips = data(1:N,8);
rakes = data(1:N,9);
quality = [4 5 5 3 3 5 5 3 3 3 4 4 5 4 5 5 4 5 5 5 3 5 3 3 4 5 3 5 5 3 3 5 3 ...
2 3 3 4 3 3 5 4 3 3 5 5 4 4 3 5 3 5 5 4 2 3 4 5 5 3 3 5 5 5 3 5]';


xs_n = xs(quality >= quality_threshold);
ys_n = ys(quality >= quality_threshold);
zs_n = zs(quality >= quality_threshold);
strikes = strikes(quality >= quality_threshold);
dips = dips(quality >= quality_threshold);
rakes = rakes(quality >= quality_threshold);
    
NN = length(strikes);
% determine the geometries of auxiliary planes
[aux_strikes, aux_dips, aux_rakes]= auxiliary_fault_plane_based_on_Kumar_and_Zhao(strikes, dips, rakes);

for i = 1:NN
if (aux_strikes(i) < 270 && aux_strikes(i) > 90)
    aux_strikes(i) =0;
    aux_dips(i) = 0;
    aux_rakes(i) = 0;
end
end

for i = 1:NN
if (strikes(i) < 270 && strikes(i) > 90)
    strikes(i) =0;
    dips(i) = 0;
    rakes(i) = 0;
end
end


for i = 1:NN
    if (strikes(i) == 0 || aux_strikes(i) == 0)
        strikes_n(i) = strikes(i)+aux_strikes(i);
        dips_n(i) = dips(i)+aux_dips(i);
        rakes_n(i) = rakes(i)+aux_rakes(i);
    elseif strikes(i) == min(strikes(i),aux_strikes(i))
        strikes_n(i) = strikes(i);
        dips_n(i) = dips(i);
        rakes_n(i) = rakes(i);
    else
        strikes_n(i) = aux_strikes(i);
        dips_n(i) = aux_dips(i);
        rakes_n(i) = aux_rakes(i);
    end
end

% strikes_n = strikes+aux_strikes;
% dips_n = dips+aux_dips;
% rakes_n = rakes+aux_rakes;


% Plotting PTI of the FMs
grid_search_results = [strikes_n' dips_n' rakes_n'];






PTIfocsphere_for_splitting_in_OADC_only_P_axis(grid_search_results);
%PTIfocsphere_for_splitting_in_OADC(grid_search_results);

% Plotting focal planes of the FMs
%figure;
%plot_strike_and_dip(grid_search_results) 

% Plot strike vs dip
figure;
plot(strikes_n,dips_n,'ro')

% Plot hypocenters
figure;plot3(xs,ys,zs,'ro')

nbins= 20;
figure;subplot(1,2,1);hist(strikes_n,nbins); subplot(1,2,2);hist(dips_n,nbins)
