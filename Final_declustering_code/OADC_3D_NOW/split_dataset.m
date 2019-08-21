



xx= xt(2,1:Nt(1));
yy= yt(2,1:Nt(1));
zz= zt(2,1:Nt(1));


outfile = 'cluster3.txt';

% write synthetic data to outfile
fid=fopen(outfile,'w');
for kk=1:Nt(1)
    
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[xx(kk) yy(kk) zz(kk)]);
    
end