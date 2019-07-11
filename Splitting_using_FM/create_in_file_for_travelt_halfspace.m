function create_in_file_for_travelt_halfspace...
    (wave, distance, depth, vel_model, halfspace)

p_incr = 0.0002; 
max_dist = round(max(distance))+500; % max distance plus 500 km. 
% This is because we may not have any ray at the maximum distance 
% based on the chosen ray parameter increment.
%
no_layers = vel_model(end,5)+1;
%
%   Now write to an ascii file
%
fid=fopen(['CSZ_travelt_',wave,'.in'],'w');
%
fprintf(fid,'%s \n','CSZ Layered Velocity Model');
fprintf(fid,'%s%s \n','CSZ_Layered_',wave);
fprintf(fid,'%s \n','   .f   .f   .f   .f   .f   .f   .f   .f');
fprintf(fid,'%s %s \n',num2str(p_incr),'  0  0');
fprintf(fid,'%s %s %s \n','0.   ',num2str(depth),'   1');
fprintf(fid,'%s %s\n','0.   ',num2str(max_dist));
fprintf(fid,'%s \n','0.');
fprintf(fid,'%s \n',num2str(no_layers)); % number of layers
%
% first layer. Atmosphere
fprintf(fid,'%s \n','  2');
fprintf(fid,'%s \n','0.001   0.001   0.001   0.0');
fprintf(fid,'%s \n','0.001   0.001   0.001   0.0');
%
for layer = 1:no_layers-1
    vp_layer = vel_model(vel_model(:,5) == layer,1);
    vs_layer = vel_model(vel_model(:,5) == layer,2);
    rho_layer = vel_model(vel_model(:,5) == layer,3);
    depth_layer = vel_model(vel_model(:,5) == layer,4);
    %
    no_el_param_in_layer = length(vp_layer);
    %
    % second layer  %'%s \n','6.00    3.50   2.7    40.0');
    fprintf(fid,'%s \n',num2str(no_el_param_in_layer));
    %
    for el_param = 1:no_el_param_in_layer    
        fprintf(fid,'%s    %s   %s     %s \n',num2str(vp_layer(el_param)),...
            num2str(vs_layer(el_param)),...
            num2str(rho_layer(el_param)),num2str(depth_layer(el_param)));
    end  
end
%
if halfspace ==1;
    if strcmp(wave,'p')
        fprintf(fid,'%s \n','  1');
        fprintf(fid,'%s \n',' 1 2');
        fprintf(fid,'%s \n',' -1 5');

    elseif strcmp(wave,'sv')    
        fprintf(fid,'%s \n','  1');
        fprintf(fid,'%s \n',' 1 2');
        fprintf(fid,'%s \n',' -1 3');

    elseif strcmp(wave,'sh')
        fprintf(fid,'%s \n','  1');
        fprintf(fid,'%s \n',' 1 2');
        fprintf(fid,'%s \n',' -1 4');

    else    
        fprintf(fid,'%s \n','Error: Specify the wave');
    end
%       
else      
    if strcmp(wave,'p')
        fprintf(fid,'%s \n','  3');
        fprintf(fid,'%s \n',' 1 2');
        fprintf(fid,'%s \n',' -1 5');
        fprintf(fid,'%s \n',' 2 2 2');
        fprintf(fid,'%s \n',' 1 8 5 ');
        fprintf(fid,'%s \n',' 2 2 2');
        fprintf(fid,'%s \n',' 1 5 5 ');
        
    elseif strcmp(wave,'sv')    
        fprintf(fid,'%s \n','  3'); % I removed the S to P wave
        fprintf(fid,'%s \n',' 1 2');
        fprintf(fid,'%s \n',' -1 3');
        fprintf(fid,'%s \n',' 2 2 2');
        fprintf(fid,'%s \n',' 1 6 3 ');
        fprintf(fid,'%s \n',' 2 2 2');
        fprintf(fid,'%s \n',' 1 3 3 ');

    elseif strcmp(wave,'sh')
        fprintf(fid,'%s \n','  3');
        fprintf(fid,'%s \n',' 1 2');
        fprintf(fid,'%s \n',' -1 4');
        fprintf(fid,'%s \n',' 2 2 2');
        fprintf(fid,'%s \n',' 1 7 4 ');
        fprintf(fid,'%s \n',' 2 2 2');
        fprintf(fid,'%s ',' 1 4 4 ');

    else    
        fprintf(fid,'%s \n','Error: Specify the wave');
    end
end


% SAMPLE
% Catchings CSV Layered Velocity Model
% CSV_layered
%    .f   .f   .f   .f   .f   .f   .f   .f
% 0.0001  0  0
% 0.   24.5   1
% 0.   500.0
% 0.
%   3
%   2
% 0.001   0.001   0.001   0.0
% 0.001   0.001   0.001   0.0
%   4
% 6.08    3.51   2.737    0.0
% 6.25    3.60   2.756    6.0
% 6.55    3.70   2.789    12.0
% 7.2     4.2     2.9     40.0
% 2
% 8.0    4.6     3.2     40.
% 8.0    4.6     3.2     200.
%  3
%  1 2
%  -1 5
%  2 2 2
%  1 8 5 
%  2 2 2
%  1 5 5
%wave = 'SV';
%depth = 25; % km
%distance = [9.424156 25.1711 35.6124 44.78913 65.8215 68.58533 331.5408 332.98 367.2564];