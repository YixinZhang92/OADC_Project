%function create_datapoints_with_FM(infile,outfile, quality_threshold)
clear all; close all; clc;

infile = 'Combined_Dataset_MT_PL.csv';
outfile = 'FM_dataset.csv';
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

xs_nn = xs_n((strikes_n ~= 0) & (dips_n ~= 0) & (rakes_n ~= 0));
ys_nn = ys_n((strikes_n ~= 0) & (dips_n ~= 0) & (rakes_n ~= 0));
zs_nn = zs_n((strikes_n ~= 0) & (dips_n ~= 0) & (rakes_n ~= 0));
strikes_nn = strikes_n((strikes_n ~= 0) & (dips_n ~= 0) & (rakes_n ~= 0));
dips_nn = dips_n((strikes_n ~= 0) & (dips_n ~= 0) & (rakes_n ~= 0));
rakes_nn = rakes_n((strikes_n ~= 0) & (dips_n ~= 0) & (rakes_n ~= 0));

% Plotting PTI of the FMs
grid_search_results = [xs_nn ys_nn zs_nn strikes_nn' dips_nn' rakes_nn'];

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:length(xs_nn)
    
    fprintf(fid,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',[xs_nn(kk) ys_nn(kk) zs_nn(kk) strikes_nn(kk) dips_nn(kk) rakes_nn(kk)]);
    
end

fclose(fid);

return