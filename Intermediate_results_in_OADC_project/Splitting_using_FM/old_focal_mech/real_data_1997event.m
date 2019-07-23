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
depth = 11.3; %11.3; %24.5; 
depth_error = 1; halfspace = 0; 

% syntax: vel_model =[vp vs rho depth layer_no]
vel_model =[6.08   3.51   2.732    0.0    1;
            6.25   3.60   2.756    6.0   1; 
            6.55   3.70   2.789    12.0   1;
            7.2    4.20   2.90    40.0   1;
            8.0    4.60   3.2    40.0001   2;
            8.0    4.60   3.2    200.0   2]; 
                
                
sta_no = [1      2       3      4       5      6     ];
dist   = [52.48  23.92   16.29  45.33   14.6   17.51 ]; 
az     = [205.35 200.21  76.28  238.16  280.4  2.11  ]; 

% Syntax: data = [station_no p_amp weight];
p_data = [1     -2094       5; % displacements are in nm
          2     2997        5;
          3     -3.178e5    5;
          4     -1.528e4    5;
          5     -2.392e4    5;
          6     294.2       5];

sh_data = [1    6165        5; %orig pos
           2    -1.196e5    5;
           3    -3.847e5    5;
           4    -3.469e4    5;
           5    1.087e5     5;
           6    -1.06e5     5];

sv_data = [1    4765        5;
           2    -2.221e4    5;
           3    8.635e4     5;
           4    6982        5;
           5    9.667e4     5;
           6    -4.002e4    5];