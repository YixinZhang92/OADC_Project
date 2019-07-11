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
vel_model =[6.55   3.70   2.789    0.0   1;
            6.55   3.70   2.789    40.0   1;
            8.0    4.60   3.2    40.0001   2;
            8.0    4.60   3.2    200.0   2]; 
                
                
sta_no = [1             2             3             4             5              6             7            ];
dist   = [6.331526e+00  1.640185e+01  1.958952e+01  1.978474e+01  3.401662e+01   4.074154e+01  4.505186e+01]; 
az     = [2.710365e+02  2.321200e+02  1.167252e+02  3.477872e+01  1.742638e+02   4.007827e+01  6.699249e+01]; 

% Syntax: data = [station_no p_amp weight];
p_data = [1     -1.185e4    5; % displacements are in nm
          2     -2326       5;
          3     2991        5;
          4     -278.6       5;
          5     1446        5;
          6     64.78       5;
          7     46.09      5];

p_data(:,2) = -1 * p_data(:,2);  

sh_data = [1    1.046e4     5;
           2    1.626e4     5;
           3    -9765       5;
           4    -235.8       5;
           5    3329        5;
           6    1137        5;
           7    -516.2      5];

sv_data = [1    7102        5;
           2    642.4       5;
           3    248.6      5;
           4    826.3       5;
           5    116.5      5;
           6    145.7       5;
           7    -165.9      5];