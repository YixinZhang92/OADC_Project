function [p_pol, sv_pol, sh_pol, vp, vs, dens] = synth_polarities...
    (dist,az,depth, halfspace, vel_model, strike, dip, rake, filename,...
     path_to_SAC, path_to_travelt, z_displ, r_displ, t_displ)
%
%   make a synthetic data set to test the focal grid search program
%
% Documentations
%
% path_to_SAC = '/usr/local/sac/bin';
% path_to_travelt = '/Users/oluwaseunfadugba/Documents/Waveform_Modeling_Now/time';
% 
[p_az, p_inc,p_eps,  ~, ~, ~] = calc_inc_from_travelt...
    (dist,az,depth, vel_model,halfspace, 'p', path_to_SAC, path_to_travelt);         
         
% [sv_az, sv_inc, ~, ~, ~] = calc_inc_from_travelt...
%     (dist,az,depth, vel_model,halfspace, 'sv', path_to_SAC, path_to_travelt);
         
[s_az, s_inc, s_eps, vp, vs, dens] = calc_inc_from_travelt...
    (dist,az,depth, vel_model,halfspace, 'sh', path_to_SAC, path_to_travelt);         
%
%  character info
ldat=length(az);
stat = repmat('stanm',ldat,1);   %eps= ones(1, ldat);
pcomp = repmat('v',ldat,1);    svcomp = repmat('r',ldat,1); 
shcomp = repmat('v',ldat,1);   pphase = repmat(' P',ldat,1); 
svphase = repmat('SV',ldat,1); shphase = repmat('SH',ldat,1);
%
% determine the moment tensor components
[m]=dismom(strike,dip,rake);
%
%  Make P wave amplitudes
%
for k=1:ldat
pamp(k)=prad(p_az(k),p_inc(k),vp,vs,dens,p_eps(k),m, z_displ);
end
%
%  Make SV wave amplitudes
%
for k=1:ldat
[svamp, svamp_r(k)] =svrad(s_az(k),s_inc(k),vp,vs,dens,s_eps(k),m, r_displ);
end
%
svamp = svamp_r; % we calculate and need only the radial component of SV wave amplitude.
%
%  Make SH wave amplitudes
%
for k=1:ldat
shamp(k)=shrad(s_az(k),s_inc(k),vp,vs,dens,s_eps(k),m, t_displ);
end
%
%   Now write to an ascii file
%
fid=fopen(filename,'w');
%
%        P wave synthetic data
%
time_p=p_inc;
for k=1:ldat
    fprintf(fid,'%-5d',[k]);
    fprintf(fid,' %5s',stat(k,1:5));
    fprintf(fid,' %5s',pcomp(k,1));
    fprintf(fid,' %15.5e %15.5e',[p_az(k) dist(k)]);
    fprintf(fid,' %5s',pphase(k,1:2));
    fprintf(fid,' %15.5e %15.5e\n',[pamp(k) time_p(k)]);
end
%
%       SV wave synthetic data
%
time_sv=s_inc;
for k=1:ldat
    fprintf(fid,'%-5d',[k+14]);
    fprintf(fid,' %5s',stat(k,1:5));
    fprintf(fid,' %5s',svcomp(k,1));
    fprintf(fid,' %15.5e %15.5e',[s_az(k) dist(k)]);
    fprintf(fid,' %5s',svphase(k,1:2));
    fprintf(fid,' %15.5e %15.5e\n',[svamp(k) time_sv(k)]);
end
%
%      SH wave synthetic data
%
time_sh = s_inc;
for k=1:ldat
    fprintf(fid,'%-5d',[k+28]);
    fprintf(fid,' %5s',stat(k,1:5));
    fprintf(fid,' %5s',shcomp(k,1));
    fprintf(fid,' %15.5e %15.5e',[s_az(k) dist(k)]);
    fprintf(fid,' %5s',shphase(k,1:2));
    fprintf(fid,' %15.5e %15.5e\n',[shamp(k) time_sh(k)]);
end

fclose(fid);

DATA_READ_FLAG=1; % flag for data read operation

%    Read a file with phase data 
fid=fopen(filename,'r');
%fid=fopen('synth_arrival_CSZ.dat','r');

%  read as characters
[B,count]=fscanf(fid,'%88c',[88 inf]);
fclose(fid);
A=B';
%  take apart the arrival
for ik=1:count
    DUMG(ik)=str2num(A(ik,1:2));
    STATG(ik,1:5)=A(ik,7:11);
    COMPG(ik,1)=A(ik,17);
    AZG(ik)=str2num(A(ik,21:33));
    DISTG(ik)=str2num(A(ik,37:49));
    PHASEG(ik,1:2)=A(ik,54:55);
    AMPG(ik)=str2num(A(ik,59:71));
    TIMEG(ik)=str2num(A(ik,75:87));
end
INCG=TIMEG;
NRECG=count;

% Setting the variables as global
% Load polarity data
    jp=0;
    jsv=0;
    jsh=0;
if DATA_READ_FLAG == 1
    k=NRECG;    % find numbers of data
    for j=1:k
        if char(PHASEG(j,1:2)) == ' P'
            jp=jp+1;
            p_pol(jp,1)=AMPG(j);%/abs(AMPG(j)); % I want the amplitudes
            p_pol(jp,2)=AZG(j);
            p_pol(jp,3)=INCG(j);
        end
        if char(PHASEG(j,1:2)) == 'SV'
            jsv=jsv+1;
            sv_pol(jsv,1)=AMPG(j);%/abs(AMPG(j));
            sv_pol(jsv,2)=AZG(j);
            sv_pol(jsv,3)=INCG(j);
        end
        if char(PHASEG(j,1:2)) == 'SH'
            jsh=jsh+1;
            sh_pol(jsh,1)=AMPG(j);%/abs(AMPG(j));
            sh_pol(jsh,2)=AZG(j);
            sh_pol(jsh,3)=INCG(j);
        end
    end
else
    p_pol=[ 0 0 0 ];
    sv_pol=p_pol;
    sh_pol=p_pol;
end
if jp == 0
    p_pol = [ 0 0 0 ];
end
if jsv == 0
    sv_pol=[ 0 0 0 ];
end
if jsh == 0
    sh_pol=[ 0 0 0 ];
end

% All amplitudes have a weight of 5. This is synthetic.
weight_P_amp = 5 * ones(length(az),1);
weight_S_amp = 5 * ones(length(az),1);

p_pol(:,4)  = weight_P_amp;
sv_pol(:,4) = weight_S_amp;
sh_pol(:,4) = weight_S_amp;

p_pol(:,5)  = (1:length(az))';
sv_pol(:,5) = (1:length(az))';
sh_pol(:,5) = (1:length(az))';

p_pol(:,6)  = (p_eps);
sv_pol(:,6) = s_eps;
sh_pol(:,6) = s_eps;

% END












% Determine the incidence angles, source parameters, and update the azimuth information. 
% Syntax: [az_now, inc_now, vp, vs, dens] = focal_error_analyses_focal_depth...
%     (x_intp,az,depth, vel_model,halfspace, wave)

% 
% clear all; close all; clc;
% %Inputs
% %  Set up azimuth and eps arrays, and determine the incidence angles
% dist = [3.5265 20 34.6410 9.3262 16.7820 28.5630 5.3590 11.5470 ...
%     14.0042 20 9.3262 20 34.6410 7.2794];
% az = [10 15. 25. 50. 80. 105. 158. 190. 210. 250. 270. 300. 320. 333.];
% eps=[ 1.  1.  1.  1.  1.   1.   1.   1.   1.   1.   1.   1.   1.   1.  ];
% 
% depth = 20; halfspace =0;
% % syntax: vel_model =[vp vs rho depth layer_no]
% vel_model =[6.08   3.51   2.732    0.0    1;
%             6.25   3.60   2.756    6.0   1; 
%             6.55   3.70   2.789    12.0   1;
%             7.2    4.20   2.90    40.0   1;
%             8.0    4.60   3.2    40.0001   2;
%             8.0    4.60   3.2    200.0   2]; 
%     
% %  Set up source orientation
% strike=30.; dip=50.; rake=55.;
% filename = 'synth_arrival.dat';
% 
% stat=['sta1 ';'sta2 ';'sta3 ';'sta4 ';'sta5 ';'sta6 ';'sta7 ';'sta8 ';'sta9 ';'sta10'; ...
%         'sta11';'sta12';'sta13';'sta14'];
% pcomp=['v';'v';'v';'v';'v';'v';'v';'v';'v';'v';'v';'v';'v';'v';];
% svcomp=['r';'r';'r';'r';'r';'r';'r';'r';'r';'r';'r';'r';'r';'r';];
% shcomp=['t';'t';'t';'t';'t';'t';'t';'t';'t';'t';'t';'t';'t';'t';];
% pphase=[' P';' P';' P';' P';' P';' P';' P';' P';' P';' P';' P';' P';' P';' P';];
% svphase=['SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV';'SV'];
% shphase=['SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH';'SH'];
% 