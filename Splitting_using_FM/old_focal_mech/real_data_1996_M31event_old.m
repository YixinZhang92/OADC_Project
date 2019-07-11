% Real data
% Inputs 
use_sv_pol = 1; % 0 -- do not use sv polarity
                % 1 -- use sv polarity
use_sv_amp = 3; % use_sv_amp = 0 -- do not use sv at all
                % use_sv_amp = 1 -- use only |sv|/sh
                % use_sv_amp = 2 -- use only |sv|/p
                % use_sv_amp = 3 -- use both |sv|/sh and |sv|/p
w_pol = 0.9; w_ratio = 0.1;
depth = 12.83; %11.3; %24.5; 
depth_error = 4; halfspace = 0; 
z_displ = 1; r_displ = 1; t_displ = 1;

% syntax: vel_model =[vp vs rho depth layer_no]
vel_model =[6.08   3.51   2.732    0.0    1;
            6.25   3.60   2.756    6.0   1; 
            6.55   3.70   2.789    12.0   1;
            7.2    4.20   2.90    40.0   1;
            8.0    4.60   3.2    40.0001   2;
            8.0    4.60   3.2    200.0   2]; 
        
        
sta_no = [1             2             3             4             5              6             7            ];
dist   = [3.401662e+01  1.958952e+01  4.505186e+01  1.640185e+01  1.978474e+01   4.074154e+01  6.331526e+00 ];
az     = [1.742638e+02  1.167252e+02  6.699249e+01  2.321200e+02  3.477872e+01   4.007827e+01  2.710365e+02 ];

% Syntax: data = [station_no p_amp weight];
p_data = [1     -3.48406e2   5;
          2     -3.54463e2   5;
          3     29.9481      5;
          4     59.32707     5;
          5     29.64097     5;
          6     -23          5;
          7     7.10383e2    5];

sv_data = [1    -66.43       5;
           2    1343.46      5;
           3    32.94        5;
           4    -395.7       5;
           5    150.1        5;
           6    60           5; 
           7    506          5];

sh_data = [1    201.4        5;
           2    -2991        5;
           3    -373.9       5;
           4    1310         5;
           5    801.9        5;
           6    217.9        5;
           7    2737         5];
