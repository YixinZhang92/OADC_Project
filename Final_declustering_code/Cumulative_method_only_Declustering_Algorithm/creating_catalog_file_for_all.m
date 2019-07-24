clear all; close all; clc;

!awk '(NR!=1) {print $4, $3, -$5}' catalog_not_in_PL17_y88_99.txt > cat_1.txt
!awk '{if ($1 <="2017") {print $6, $5, -$7}}' catalog_not_in_PL17_y00_17.txt  >> cat_1.txt

infile = 'cat_1.txt';

fid=fopen(infile,'r');
[data1,~]=fscanf(fid,'%g %g %g',[3,inf]);
fclose(fid);

data1=data1'; 
    
data1n(:,1) = (data1(:,1)+70.6) *75.778; 
data1n(:,2) = (data1(:,2)-47.2) *111.1743; 
data1n(:,3) = data1(:,3);


data1n = data1n(data1n(:,1)>-10 & data1n(:,2)>-10 & data1n(:,1)< 80 & data1n(:,2)< 100 ,:);




!awk '$0 !~ /Station/ {print $1, $2, $3}' CSZ_hypos.txt > cat_2.txt

infile2 = 'cat_2.txt';

fid=fopen(infile2,'r');
[data2,~]=fscanf(fid,'%g %g %g',[3,inf]);
fclose(fid);

data2=data2'; 


data = [data1n ; data2]; [nhypos,~]=size(data);





% xmin = min(data(:,1));
% ymin = min(data(:,2));
% 
% 
% data(:,1) = data(:,1) - xmin + 10;
% data(:,2) = data(:,2) - ymin + 10;
% 

outfile = 'All_CSZ_hypos.txt';

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:nhypos
    
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[data(kk,1) data(kk,2) data(kk,3)]);
    
end


figure; plot3(data(:,1),data(:,2),data(:,3),'ro'); grid MINOR


% plot3(data1n(:,1),data1n(:,2),data1n(:,3),'bo'); hold on; 
% plot3(data2(:,1),data2(:,2),data2(:,3),'ro');
