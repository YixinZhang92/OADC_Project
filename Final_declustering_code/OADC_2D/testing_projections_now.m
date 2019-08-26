clear all; close all; clc;
global xs ys zs N

infile = 'Simul.1_hypos.txt';
%infile = 'All_CSZ_hypos.txt'; 
%infile = 'Simul.3_ALL_hypos_hypos.txt'
simul_tag = 'Simul1_hypos.';
read_catalog(infile,simul_tag); %view(-az,el);
X=[xs' ys' zs'];
        
%az = 30;
%el = 30;

for az=30%0:5:360
    for el=60
     
        U = [cosd(az), -sind(az), 0];
        V = [-sind(az)*sind(el), -cosd(az)*sind(el), -cosd(el)];
        n = [sind(az)*cosd(el), cosd(az)*cosd(el), -sind(el)];

        R = [U' V' n'];
        projX=R\X';
        pxs = projX(1,:);
        pys = -projX(2,:);

        Fig1 = figure('Name','Projected hypocenters','Position', get(0, 'Screensize'));
        subplot(1,2,1)
        plot3(xs,ys,zs,'o');axis equal; view(-az, el);grid MINOR

        subplot(1,2,2)
        plot(pxs,pys,'o');grid MINOR; shg

    end
end

outfile = [simul_tag 'proj.' num2str(az) '.' num2str(el) '.txt'];

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:N
    
    fprintf(fid,'%12.5f %12.5f\n',[pxs(kk) pys(kk)]);
    
end

