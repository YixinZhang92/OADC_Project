% Real data
z_displ = 0; r_displ = 1; t_displ = 1;

% Inputs 
use_sv_pol = 1; % 0 -- do not use sv polarity
                % 1 -- use sv polarity
use_sv_amp = 3; % use_sv_amp = 0 -- do not use sv at all
                % use_sv_amp = 1 -- use only |sv|/sh
                % use_sv_amp = 2 -- use only |sv|/p
                % use_sv_amp = 3 -- use both |sv|/sh and |sv|/p
w_pol = 0.8; w_ratio = 0.2;
depth = 12.83; %11.3; %24.5; 
depth_error = 4; halfspace = 0; 

% syntax: vel_model =[vp vs rho depth layer_no]
vel_model =[6.08   3.51   2.732    0.0    1;
            6.25   3.60   2.756    6.0   1; 
            6.55   3.70   2.789    12.0   1;
            7.2    4.20   2.90    40.0   1;
            8.0    4.60   3.2    40.0001   2;
            8.0    4.60   3.2    200.0   2]; 
        
        
sta_no = [1             2             3             4             5              6             7            ];
dist   = [6.331526e+00  1.640185e+01  1.958952e+01  1.978474e+01  3.401662e+01   4.074154e+01  4.505186e+01 ];
az     = [2.710365e+02  2.321200e+02  1.167252e+02  3.477872e+01  1.742638e+02   4.007827e+01  6.699249e+01 ];

% Syntax: data = [station_no p_amp weight];
p_data = [1     -11.5536     5;
          2     -0.49773     5;
          3     2.64871      5;
          4     -0.57079     5;
          5     2.01622      5;
          6     0.265527     5;
          7     0.0252386    5];
      
p_data(:,2) = -1 * p_data(:,2);      
       
% sv_data = [1    0.27243      5; % from polarization angle
%            2    5.676        5;
%            3    0.41315      5;
%            4    -2.31616     5;
%            5    -0.03226     5;
%            6    0.955088     5;
%            7    0.08938      5];

sv_data = [1    -11.57       5; % by visualization
           2    -3.413       5;
           3    -5.448       5;
           4    -3.634       5;
           5    -0.5456      5;
           6    0.8538       5;
           7    0.1648       5];
% sv_data(:,2) = -1 * sv_data(:,2);      
      
sh_data = [1    -12.23       5; 
           2    26.16        5;
           3    -19.96       5;
           4    1.313        5;
           5    4.638        5;
           6    4.233        5; 
           7    -0.3843      5];