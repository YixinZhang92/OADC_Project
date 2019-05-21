close all; clear all; clc;
tic

%infile = 'testdata.txt';
infile = 'declustered_collapsed_hypo_0.6.txt';
%infile= 'result_declustered_collapsed_hypo_0.8.txt';
%infile= 'result_collapsed_hypo_0.8.txt';
%infile= 'result_declustered_hypo_0.8.txt';
global xs ys zs 

%***************** Read Catalog of Hypocenters ****************************
read_catalog(infile);

% Random hypocenters have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike

% rotate into rake direction
nhypos = length(xs);
strike_incr = 5; dip_incr= 5;

%The solution is not sensitive to angles greater than 180. 181 deg is the same as 1 deg
strikes = 0:strike_incr:90;%175;%[45 135];%pos is anticlockwise. 
dips = 0:dip_incr:90;%[90 45];

width= 1; % width of the depth slider when determining the number of EQs in each block.
min_eqs_for_a_cluster = 20;
mult_incr = 1;

% Initializing array %analy_old = analy;
%7 - no of eqs, 8 - dist at max eqs
analy = [xs' ys' zs' zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1) ...
    zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1)  zeros(nhypos,1)]; 

unique_value = 0;

for j = 1: length(strikes)
    for jj = 1:length(dips)
        
        % Getting a unique value
        unique_value = unique_value + 1; 
        
        % ------------------------------------------------------------------------
        % Analyzing for each combination of strike and dip. We do not need the
        % rake.
        con=pi/180.;
        strike=strikes(j).*con;
        dip=dips(jj).*con;

        % extracting the xs, ys and zs because their relative location will change
        % during the analysis.
        xs_now = analy(:,1)';
        ys_now = analy(:,2)';
        zs_now = analy(:,3)';

        xs_now = xs_now - mean(xs_now); % removing the mean for the rotation to be about origins. 
        ys_now = ys_now - mean(ys_now);
        zs_now = zs_now - mean(zs_now);

        R=[xs_now ; ys_now ;zs_now];

        % rotate into strike direction
        %Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
        Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];

        Rstrike=Dstrike*R;

        rxp(1:nhypos) = Rstrike(1,1:nhypos);
        ryp(1:nhypos) = Rstrike(2,1:nhypos);
        rzp(1:nhypos) = Rstrike(3,1:nhypos);

        % rotate into dip direction
        Rstrike(1,1:nhypos) = Rstrike(1,1:nhypos) - mean(Rstrike(1,1:nhypos));
        Rstrike(2,1:nhypos) = Rstrike(2,1:nhypos) - mean(Rstrike(2,1:nhypos));
        Rstrike(3,1:nhypos) = Rstrike(3,1:nhypos) - mean(Rstrike(3,1:nhypos));

        % Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
        %Ddip=[ cos(dip) 0 sin(dip); 0 1  0 ; -sin(dip) 0 cos(dip)];
        Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];

        Rdip=Ddip*Rstrike;

        Rdip(1,1:nhypos) = Rdip(1,1:nhypos) + mean(Rstrike(1,1:nhypos));
        Rdip(2,1:nhypos) = Rdip(2,1:nhypos) + mean(Rstrike(2,1:nhypos));
        Rdip(3,1:nhypos) = Rdip(3,1:nhypos) + mean(Rstrike(3,1:nhypos));

        rxp(1:nhypos) = Rdip(1,1:nhypos);
        ryp(1:nhypos) = Rdip(2,1:nhypos);
        rzp(1:nhypos) = Rdip(3,1:nhypos);

        % Take the depth to have a minimu of zero. Thi does not overright the
        % original dataset.
        rzp = rzp-min(rzp);
        
        % Concatenate these result with the original dataet to kep track of each
        % point. Sort the data based on the new depths.
        analy(:,4) = rxp';
        analy(:,5) = ryp';
        analy(:,6) = rzp';

        analy = sortrows(analy,6);

        % Stepwisely move through the data in depth and the number of EQs in each
        % block. we need to specify the width of the block. 
        st_array = 0:max(rzp+1); 

        for i = 1:length(st_array)
            st = st_array(i);
            index = find(analy(:,6)>=st & analy(:,6)<(st+width));

            neqs_in_window = length(index);

            for ii = 1: neqs_in_window
                if (analy(index(ii),7) < neqs_in_window) && ...
                        (neqs_in_window > min_eqs_for_a_cluster)
                    analy(index(ii),7) = neqs_in_window;
                    analy(index(ii),8) = st;
                    analy(index(ii),9) = strike/con;
                    analy(index(ii),10) = dip/con;
                    analy(index(ii),11) = unique_value;
                end
            end    
        end
    end
end

analy = sortrows(analy,11);

[clus,ia,ic] = unique(analy(:,11));
a_counts = accumarray(ic,1);
value_counts = [clus, a_counts];
value_counts = sortrows(value_counts,2,'descend');

for i= 1:length(value_counts(:,1))
    dist_index = analy(analy(:,11)== value_counts(i,1),8); value_counts(i,3) = dist_index(1);
    strike_a = analy(analy(:,11)== value_counts(i,1),9); value_counts(i,4) = strike_a(1);
    dip_a = analy(analy(:,11)== value_counts(i,1),10); value_counts(i,5) = dip_a(1);

end

value_counts = round(value_counts);

value_counts(:,6) = 0;

for i= 1:length(value_counts(:,1))
    
    if value_counts(i,6) == 0
        
        value_counts((abs(value_counts(:,3) - value_counts(i,3))<=1) & ...
            (abs(value_counts(:,4) - value_counts(i,4))<=mult_incr*strike_incr) & ...
            (abs(value_counts(:,5) - value_counts(i,5))<=mult_incr*dip_incr)& ...
            (value_counts(:,6) == 0) ,6) = i;

    end
end


uni = unique(value_counts(:,6));

for i=1:length(uni)
    uniq_index = value_counts(value_counts(:,6) == uni(i),1);

    for j = 1:length(uniq_index)
        % I need to replace the uniq_index last.
        analy(analy(:,11) == uniq_index(j),9) = value_counts(uni(i),4); % strike
        analy(analy(:,11) == uniq_index(j),10) = value_counts(uni(i),5); % dip
        analy(analy(:,11) == uniq_index(j),8) = value_counts(uni(i),3); % dist
        analy(analy(:,11) == uniq_index(j),7) = sum(value_counts(value_counts(:,6) == uni(i),2)); % no_of_eqs----
        
        analy(analy(:,11) == uniq_index(j),11) = value_counts(uni(i),1); % strike
        
    end
    
    
end


analy = sortrows(analy,11);

[clus,ia,ic] = unique(analy(:,11));
a_counts = accumarray(ic,1);
value_counts = [clus, a_counts];
value_counts = sortrows(value_counts,2,'descend');

cluster = value_counts(:,1);

% ------------------------------------------------------------------------
% Creating figures
% ------------------------------------------------------------------------
fig = figure;
ax1 = subplot(1,2,1);
plot3(xs,ys,zs,'o');
axis equal; title('Input Hypocenters');
xlabel('X km'); ylabel('Y km'); zlabel('Z km'); grid on

ax2 = subplot(1,2,2); i = 0;
for ncluster = 1:length(cluster)
    
    if value_counts(ncluster,2) >= min_eqs_for_a_cluster
        xsc = analy(analy(:,11) == cluster(ncluster),1);
        ysc = analy(analy(:,11) == cluster(ncluster),2);
        zsc = analy(analy(:,11) == cluster(ncluster),3);

        plot3(xsc,ysc,zsc,'o'); hold on;
        
        i = i+1;
    
    else
        xsc = analy(analy(:,11) == cluster(ncluster),1);
        ysc = analy(analy(:,11) == cluster(ncluster),2);
        zsc = analy(analy(:,11) == cluster(ncluster),3);

        plot3(xsc,ysc,zsc,'ko'); hold on;
    
    end
end

axis equal; title('Clustered Data'); 
xlabel('X km'); ylabel('Y km'); zlabel('Z km'); grid on

hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); 
rotate3d on

% OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)

toc 

%*************************** END ******************************************