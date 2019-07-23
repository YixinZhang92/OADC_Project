function plot_Swave_pola_foceqarea(svamp,shamp, az, inc, rot_angle_input,arrowcolor)
% Plotting the S wave polarization vector on the P wave readiation pattern.
% This is purely analytical. The displacements are calculated from the
% given fault plane solutions.

con=pi/180.; radius=1.0;
arrow_length = sqrt(svamp.^2+shamp.^2);
max_len = max(max(arrow_length)); 

for m=1:length(az);
        % samplying azimuth and incidence angle, extracting teh
        % corresponding SV and SH amplitudes.
        azn = az(m); incn = inc(m);
        svampn = svamp(m); shampn = shamp(m);
        %
        r=sqrt(2)*radius*sin(incn*con/2);
        ycnod=r.*cos(azn*con);
        xcnod=r.*sin(azn*con);
        %
        % max length to be 0.2
        scaled_length = (sqrt(svampn.^2+shampn.^2)/max_len)*0.2;
        %
        rot_angle = rot_angle_input(m);
        %
        % create arrow using necessary parameters and place it on the focal
        % sphere.
        create_arrow_2(xcnod,ycnod,scaled_length,rot_angle,arrowcolor);
end;
%
title('S-wave Polarization vector plotted on P-wave radiation pattern');

















%         % how to project the angles on the lower hemisphere
%         plot([p0_bar(1) p1_bar(1)],[p0_bar(2) p1_bar(2)],'k'); hold on;
%         plot([p1_bar(1) ah_1_bar(1)],[p1_bar(2) ah_1_bar(2)],'k'); hold on;
%         plot([p1_bar(1) ah_2_bar(1)],[p1_bar(2) ah_2_bar(2)],'k'); hold on;