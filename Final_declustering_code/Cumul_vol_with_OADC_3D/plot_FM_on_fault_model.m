function plot_FM_on_fault_model(infile_FM,quality_threshold,use_mag_size,fm_size_or_factor)
global xs ys zs

fid=fopen(infile_FM,'r');

[data,count]=fscanf(fid,'%g %g %g',[11,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xsf=data(1:N,1);
ysf=data(1:N,2);
zsf=data(1:N,3);
mags = data(1:N,6);
strikes = data(1:N,7);
dips = data(1:N,8);
rakes = data(1:N,9);
quality = [4 5 5 3 3 5 5 3 3 3 4 4 5 4 5 5 4 5 5 5 3 5 3 3 4 5 3 5 5 3 3 5 3 ...
2 3 3 4 3 3 5 4 3 3 5 5 4 4 3 5 3 5 5 4 2 3 4 5 5 3 3 5 5 5 3 5]';

% 
xs_n = xsf(quality >= quality_threshold);
ys_n = ysf(quality >= quality_threshold);
zs_n = -zsf(quality >= quality_threshold);
mags = mags(quality >= quality_threshold);
strikes = strikes(quality >= quality_threshold);
dips = dips(quality >= quality_threshold);
rakes = rakes(quality >= quality_threshold);


% Find the closet FM to the earthquake in the thick cluster
kk= 0; dst = zeros(1,length(xs_n));
for i = 1:length(xs) % per hypocenter

    for m=1:length(xs_n)  %per FM
        dst(m)= sqrt((xs_n(m)-xs(i))^2 + (ys_n(m)-ys(i))^2 + (zs_n(m)-zs(i))^2);           
    end

    %  find the closest fault plane
    [mindist, index] = min(dst);

    if mindist < 1 %dist2FM_threshold 
        kk = kk + 1;

        xs_fm(kk) = xs_n(index);
        ys_fm(kk) = ys_n(index);
        zs_fm(kk) = zs_n(index);
        mags_fm(kk) = mags(index);

        strikes_fm(kk) = strikes(index);
        dips_fm(kk) = dips(index);
        rakes_fm(kk) = rakes(index);        
    end
end

fm = [strikes_fm' dips_fm' rakes_fm'];
centerX = xs_fm;
centerY = ys_fm;
centerZ = zs_fm;

if use_mag_size ==0
    diam = ones(1,length(xs_fm))*fm_size_or_factor;
else
    diam = ones(1,length(xs_fm)).*mags_fm*fm_size_or_factor;
end

ta = 0;
color ='k';
beachball_CERI_code_3D(fm, centerX, centerY, centerZ, diam, ta, color)
