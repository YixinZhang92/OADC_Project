clear all; close all; clc;

strikes = [80 80 78 78 78 78 76 106 100 80 80 80 102 216 80 76 78 78 80 78 30 120];
dips = [30 30 30 30 30 30 30 32 28 30 30 30 26 76 28 36 30 30 32 30 22 32];
rakes = [116 116 114 114 116 114 114 154 150 116 116 116 150 68 110 126 114 114 118 114 88 144];

N = length(strikes); az=linspace(0,360,1000); inc=linspace(0,90,300); P_contour = zeros(300,1000);

for i = 1:N
%   calculate the moment tensor for a point dislocation
m=dismom(strikes(i),dips(i),rakes(i));
%
%   calculate the P, T, and I vectors for this moment tensor
svec=stressvec(m);

ind = find(svec(:,1) < 0); svec(ind,1) = 360 + svec(ind,1);

P_coord(i,:) = round(svec(1,:),1);
T_coord(i,:) = round(svec(2,:),1);    
I_coord(i,:) = round(svec(3,:),1);    
  
v = abs(az - P_coord(i,1)); az_ind = find(v == min(v));
vi = abs(inc - P_coord(i,2)); inc_ind = find(vi == min(vi)); %P_ind(i,:) = [az_ind inc_ind];
P_contour(inc_ind, az_ind) = P_contour(inc_ind, az_ind) + 1;

end
% 
%
N = max(max(P_contour));

c_all = zeros(2*N,5000);
icontour_all = zeros(N,1);
numstart_all = zeros(N,5);
numfin_all = zeros(N,5);
%
for no_of_FM = 1%:N
        %  find the zero amplitude contour, if it exists
        c=contour(az,inc,P_contour,[1 1]);
        %
        %  find all of the contours
        [ix,iy]=size(c);
        icontour=0;
        num=0;
        numstart=0;
        numfin=0;
        if iy > 0
        num1=c(2,1);
        if (num1 > 1)
            icontour=1;
            num(1)=num1;
            numstart(1)=2;
            numfin(1)=numstart(1)+num1-1;
        else
            icontour=0;
        end
        size(numfin);
        istop=0;
        ik=0;
        if (icontour > 0)
            test=numfin(1);
            while (test < iy)
                ik=ik+1;
                icontour=icontour+1;
                num(ik+1)=c(2,numfin(ik)+1);
                numstart(ik+1)=numfin(ik)+2;
                numfin(ik+1)=numstart(ik+1)+num(ik+1)-1;
                test=numfin(ik+1);
            end
        end
        end
 
        [nsr, nsc] = size(numstart);
        [nfr, nfc] = size(numfin);
        icontour_all(no_of_FM,:) = icontour;
        numstart_all(no_of_FM,1:nsc) = numstart;
        numfin_all(no_of_FM,1:nfc) = numfin;
        
        [mmm,nnn]=size(c);
        a = 2*no_of_FM - 1;        b = 2*no_of_FM;
    
        c_all(a:b,1:nnn) = c;
        
end

data_flag = 1;
line_color = 'r';
%p_pol = [p_amp  p_az        p_inc       weight_p_amp  p_index p_eps];
p_pol = [12930   271.0365    27.461801831175     5    1    -1;
         2704    232.12      55.0736556719227    5    2    -1;
         -3135   116.7252    60.4914228157975    5    3    -1;
         33.1    34.77872    60.7884642041669    5    4    -1;
         -1715   174.2638    75.762800373442     5    5    -1;
         -127.5  40.07827    80.1389504705361    5    6    -1;
         136.9   66.99249    82.416314722158     5    7    -1];


foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,p_pol, data_flag,line_color)













% 
% 
% %
% c_all = zeros(2*N,5000);
% icontour_all = zeros(N,1);
% numstart_all = zeros(N,5);
% numfin_all = zeros(N,5);
% %
% for no_of_FM = 1:N
%         eps=+1.0;
%         %   calculate the moment tensor for a point dislocation
%         m=dismom(strike(no_of_FM),dip(no_of_FM),rake(no_of_FM));
%         %
%         %  Compute P wave radiation pattern
%         az=linspace(0,360,1000);
%         inc=linspace(0,90,300);
%         pamp=prad(az,inc,vp,vs,dens,eps,m, z_displ);
%         
%         save pamp.mat pamp
%         %
%         %  find the zero amplitude contor, if it exists
%         c=contourc(az,inc,pamp,[0 0]);
%         %
%         %  find all of the contours
%         [ix,iy]=size(c);
%         icontour=0;
%         num=0;
%         numstart=0;
%         numfin=0;
%         if iy > 0;
%         num1=c(2,1);
%         if (num1 > 1);
%             icontour=1;
%             num(1)=num1;
%             numstart(1)=2;
%             numfin(1)=numstart(1)+num1-1;
%         else;
%             icontour=0;
%         end;
%         size(numfin);
%         istop=0;
%         ik=0;
%         if (icontour > 0);
%             test=numfin(1);
%             while (test < iy);  
%                 ik=ik+1;
%                 icontour=icontour+1;
%                 num(ik+1)=c(2,numfin(ik)+1);
%                 numstart(ik+1)=numfin(ik)+2;
%                 numfin(ik+1)=numstart(ik+1)+num(ik+1)-1;
%                 test=numfin(ik+1);
%              end;
%          end;
%         end;
%  
%         [nsr, nsc] = size(numstart);
%         [nfr, nfc] = size(numfin);
%         icontour_all(no_of_FM,:) = icontour;
%         numstart_all(no_of_FM,1:nsc) = numstart;
%         numfin_all(no_of_FM,1:nfc) = numfin;
%         
%         [mmm,nnn]=size(c);
%         a = 2*no_of_FM - 1;        b = 2*no_of_FM;
%     
%         c_all(a:b,1:nnn) = c;
%         
% end
