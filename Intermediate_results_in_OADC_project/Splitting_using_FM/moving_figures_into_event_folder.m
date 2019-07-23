event_no = 120; 
MT_FM = 1;

if MT_FM == 1
    !mv fig_P_pol.png fig_P_pol_MTFM.png
    !mv fig_SV_pol.png fig_SV_pol_MTFM.png
    !mv fig_SH_pol.png fig_SH_pol_MTFM.png
    
    eval(sprintf('%s%s%s%s%s%s', '!mv fig_*MTFM.png ',...
        'real_data_focal2_CSZ_no_SVpol/event',num2str(event_no),'_FM;'));

else     
    eval(sprintf('%s%s%s%s%s%s','mkdir real_data_focal2_CSZ_no_SVpol/event',...
        num2str(event_no),'_FM; !mv fig_* write_mechanisms.txt ',...
        'real_data_focal2_CSZ_no_SVpol/event',num2str(event_no),'_FM;'));
end

% eval(sprintf('%s%s%s%s%s%s','mkdir real_data_focal2_CSZ/event',...
%     num2str(event_no),'_FM; !mv fig_* write_mechanisms.txt ',...
%     'real_data_focal2_CSZ/event',num2str(event_no),'_FM;'));