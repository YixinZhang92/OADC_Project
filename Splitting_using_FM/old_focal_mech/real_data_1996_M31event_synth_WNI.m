% Real data
z_displ = 1; r_displ = 1; t_displ = 1;

% Inputs 
use_sv_pol = 1; % 0 -- do not use sv polarity
                % 1 -- use sv polarity
use_sv_amp = 3; % use_sv_amp = 0 -- do not use sv at all
                % use_sv_amp = 1 -- use only |sv|/sh
                % use_sv_amp = 2 -- use only |sv|/p
                % use_sv_amp = 3 -- use both |sv|/sh and |sv|/p
w_pol = 0.7; w_ratio = 0.3;
depth = 12.83; %11.3; %24.5; 
depth_error = 1; halfspace = 0; 

% syntax: vel_model =[vp vs rho depth layer_no]
vel_model =[6.08   3.51   2.732    0.0    1;
            6.25   3.60   2.756    6.0   1; 
            6.55   3.70   2.789    12.0   1;
            7.2    4.20   2.90    40.0   1;
            8.0    4.60   3.2    40.0001   2;
            8.0    4.60   3.2    200.0   2]; 
                
                
sta_no = [1             2             3             4             5              6             7            ];
dist   = [6.331526e+00  1.640185e+01  1.958952e+01  1.978474e+01  3.401662e+01   4.074154e+01  4.505186e+01]; 
az     = [2.710365e+02  2.321200e+02  1.167252e+02  3.477872e+01  1.742638e+02   4.007827e+01  6.699249e+01]; 

% Syntax: data = [station_no p_amp weight];
p_data = [1     -1.293e4    5; % displacements are in nm
          2     -2704       5;
          3     3135        5;
          4     -33.1       5;
          5     1715        5;
          6     127.5       5;
          7     -136.9      5];

p_data(:,2) = -1 * p_data(:,2);  

sh_data = [1    1.142e4     5;
           2    1.764e4     5;
           3    -1.026e4    5;
           4    888.9       5;
           5    3378        5;
           6    1963        5;
           7    -105.7      5];

sv_data = [1    6272        5;
           2    -1857       5;
           3    -403.8      5;
           4    734.7       5;
           5    -142.8      5;
           6    71.83       5;
           7    -294.7      5];