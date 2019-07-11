% Real data
z_displ = 1; r_displ = 1; t_displ = 1;

% Inputs 
use_sv_pol = 1; % 0 -- do not use sv polarity
                % 1 -- use sv polarity
use_sv_amp = 3; % use_sv_amp = 0 -- do not use sv at all
                % use_sv_amp = 1 -- use only |sv|/sh
                % use_sv_amp = 2 -- use only |sv|/p
                % use_sv_amp = 3 -- use both |sv|/sh and |sv|/p
w_pol = 0.8; w_ratio = 0.2;
depth = 16.0; %11.3; %24.5; 
depth_error = 1; halfspace = 0; 

% syntax: vel_model =[vp vs rho depth layer_no]
vel_model =[6.08   3.51   2.732    0.0    1;
            6.25   3.60   2.756    6.0   1; 
            6.55   3.70   2.789    12.0   1;
            7.2    4.20   2.90    40.0   1;
            8.0    4.60   3.2    40.0001   2;
            8.0    4.60   3.2    200.0   2]; 
                
                
sta_no = [1      2       3      4       5      6     ];
dist   = [64.01  38.14   29.42  47.88   13.33  11.96 ]; 
az     = [190.07 175.53  113.36 215.02  193.8  80.66 ]; 

% Syntax: data = [station_no p_amp weight];
p_data = [1     -823.5      5; % displacements are in nm
          2     286.6       5;
          3     -5522       5;
          4     -112.6      5;
          5     1.319e4     5;
          6     2.106e4     5];

sh_data = [1    4321        5; %orig pos
           2    -4270       5;
           3    25390       5;
           4    -3615       5;
           5    -3.174e4    5;
           6    4.532e4     5];

sv_data = [1    -1647       5;
           2    -8310       5;
           3    -3.07e4     5;
           4    3973        5;
           5    -6.357e4    5;
           6    -1.695e5    5];