% Real data for 1997        
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
        
             
sta_no = [1        2        3        4         5         6         7        8        9      10       11];
dist   = [5.25E+01 2.40E+01 1.66E+01 4.51E+01 1.42E+01 1.73E+01  9.33E+02 8.22E+02 8.25E+02 8.52E+03 8.52E+03];
az     = [2.05E+02 1.99E+02 7.73E+01 2.38E+02 2.80E+02 3.38E+00 1.05E+02 2.43E+02 1.39E+01 1.05E+02 1.05E+02];
            
 % Syntax: data = [station_no p_amp weight];
        p_data = [1     -371.50105      5;
                  2     1.46418e3       5;
                  3     -4.77337e4    	5;
                  4     -2.96629e3      5;
                  5     -4.369166e3     5;
                  6     477.8573        5;
                  7     12.585          5;
                  8     2.10956         5;
                  %9     18.97857        5;
                  10    -5.8315         5;
                  11    -2.19           5];

        sv_data = [1    -1.40862e3      5; % by Spol angle and observation. Stretching and compression in time was observed
                   2    -5364.62        5;
                   3    -70248.2        5;
                   4    -4.78104e3      5;
                   5    60034           5;
                   6    -2575.313       5];

        sh_data = [1    -1.8864e3       5;
                   2    -1.23305e4   	5;
                   3    -6.81192e4  	5;
                   4    -1.12957e4  	5;
                   5    107520.2        5;
                   6    -1.0681e4    	5];
               
               
               
% END               
               
               
               
               
               
               
               
               
               
               
               
               
               
               %         % Syntax: data = [station_no p_amp weight];
%         p_data = [1     -6.81895e2      5;
%                   2     1.43831e3       5;
%                   3     -4.82908e4    	5;
%                   4     -2.96568e3      5;
%                   5     -4.04418e3      5;
%                   6     5.37341e1      5              ;
%                   7     10              2;
% %                   8     2.08812         1;
% %                   9     -1.23303        1;
%                   10    -10             2;
%                   11    -5              2];
% 
% %               ;
% %                   7     10              2;
% %                   8     2.08812         2;
% %                   9     -1.23303        2;
% %                   10    -10             2;
% %                   11    -5              2
%                   
% %         sv_data = [1    1.73134e3         5; % by reading seismogram
% %                    2    3093.9945         5;
% %                    3    -1.19885e5        5;
% %                    4    -4.87407e3        5;
% %                    5    2.49798e5         5;
% %                    6    -1.12332e4        5];
%                
%         sv_data = [1    704.706         5; % by Spol angle
%                    2    3093.9945       5;
%                    3    22412.9817      5;
%                    4    -1753.7921      5;
%                    5    201515.9723     5;
%                    6    -2575.313       5];
% 
%         sh_data = [1    -1.8864e3       5;
%                    2    -1.08961e4   	5;
%                    3    -6.65931e4  	5;
%                    4    -1.12957e4  	5;
%                    5    1.5769e5        5;
%                    6    -1.12332e4    	5];
