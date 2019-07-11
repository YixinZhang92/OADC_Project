% This parameter.in file is an input for focal.m program.
%
% After editting this file, run focal.m on the MATLAB terminal using "focal"

path_to_SAC = '/usr/local/sac/bin';
path_to_travelt = '/Users/oluwaseunfadugba/Documents/Waveform_Modeling_Now/time';
% -------------------------------------------------------------------------

%  vel_model = [6.0    3.5   2.7    0.0    1; ... % Halfspace model
%               6.0    3.5   2.7    40.0   1];
%
synthetic = 1; 
%
if synthetic == 1
        % Generate and using synthetic dataset         
        % Set up source orientation
        strike= 36; dip=45.; rake=90.;
        
        filename = 'synth_arrival.dat';
        %  Set up azimuth and eps arrays, and determine the incidence angles
%          dist = [3.5265 20 34.6410 9.3262 16.7820 28.5630 5.3590 11.5470 ...
%                  14.0042 20 93262 2000 346410 72794];
%          az = [10 15. 25. 50. 80. 105. 158. 190. 210. 250. 270. 300. 320. 333.];

         dist = [3.5265 20 3400.6410 6000.3262 160.7820 280.5630 ... 
             4000.3590 110.5470 140.0042 200 932 2000  ...
                500.3590 1001.5470 1400.0042 200 8362 2000 3410 724];
         az = [10 15. 25. 50. 80. 105. 100 120 122 130 125 130 158. 190. 210. 250. 270. 300. 320. 333.];
%         dist = [6.331526 16.40185 19.58952 19.78474 34.01662 40.74154 45.05186];
%         az       = [271.0365 232.120  116.7252 34.77872 174.2638 40.07827 66.99249];
        
        z_displ = 1; r_displ = 1; t_displ = 1;
        
        % Inputs 
        use_sv_pol = 1; % 0 -- do not use sv polarity
                        % 1 -- use sv polarity
        use_sv_amp = 3; % use_sv_amp = 0 -- do not use sv at all
                        % use_sv_amp = 1 -- use only |sv|/sh
                        % use_sv_amp = 2 -- use only |sv|/p
                        % use_sv_amp = 3 -- use both |sv|/sh and |sv|/p
        w_pol = 0.7; w_ratio = 0.3;
        depth =  12.83; %11.3; %24.5; 
        depth_error = 1; halfspace = 0; 

        % syntax: vel_model =[vp vs rho depth layer_no]
        vel_model =[6.08   3.51   2.732    0.0    1;
                    6.25   3.60   2.756    6.0   1; 
                    6.55   3.70   2.789    12.0   1;
                    7.2    4.20   2.90    40.0   1;
                    8.0    4.60   3.2    40.0001   2;
                    8.0    4.60   3.2    200.0   2]; 
        
        % Note that since we only calculate or measure the radial component
        % of SV-wave amplitude, there will be a discrepancy between the
        % 'observed' and calculated S-wave polarization vectors. To remove this 
        % differences, we need to request the user to give the vertical 
        % component amplitude of the Swave as well. 
        
else
        % The inputs for the real data are supplied by the user via a
        % MATLAB function. The inputs are station number, epicentral distance
        % to the stations, azimuth of the stations, p_data, sv_data and sh_data.
        % E.g.
        %         sta_no = [1        2       ];
        %         dist   = [9.424156 25.1711 ];
        %         az     = [23.3895  126.803 ];
        %
        % Syntax: data = [station_no p_amp weight];
        %         p_data = [1     6.00804e3   5;
        %                   2     6.02727e2   5;

        
        %run real_data_focal2_CSZ/real_data_event1()
        run real_data_focal2_CSZ_no_SVpol/real_data_event120()
end


%% Run this code afterward to move figures to the right folder
% i.e. mkdir real_data_focal2_CSZ/event17_FM; !mv fig_* write_mechanisms.txt 
% real_data_focal2_CSZ/event17_FM);

% event_no = 17; 
% eval(sprintf('%s%s%s%s%s%s','mkdir real_data_focal2_CSZ/event',...
%     num2str(event_no),'_FM; !mv fig_* write_mechanisms.txt ',...
%     'real_data_focal2_CSZ/event',num2str(event_no),'_FM;'));











%         
%         real_data_event1()
%         real_data_event5()
%         real_data_event7()
%         
%         
%         
%         %real_data_1997event()
%         %real_data_1996_M31event()
%         %real_data_2013_192event()
%         %real_data_1996_M31event_synth_WNI()
%         %real_data_1996_M31event_synth_HALFSPACE_WNI()
%         
%         %real_data_event32()
%         %real_data_1996_M31event_now()
%         %real_data_wavenumber_integration_systh()