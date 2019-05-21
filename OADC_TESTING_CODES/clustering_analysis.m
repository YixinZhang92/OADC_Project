function [n0,analy, value_counts] = clustering_analysis...
    (xs, ys, zs, strike_incr, dip_incr, width, mult_incr, min_eqs_for_a_cluster)



% Inputs: strike_incr, dip_incr - strike and dip increments
%         width (in km) width of the depth slider when determining the number of EQs in each block.
%         min_eqs_for_a_cluster - min number of earthquakes allowed for a cluster
%         mult_inc - =1, multiple of strike and dip increments used when
%         combining the datasets with similar strike and dip
%
% Outs:  n0 - no of clusters with eqs >= min_eqs_for_a_cluster
%        analy = [x y z neqs st strike dip unique_value]
%
% strike_incr = 5; dip_incr= 10;
% width= 1; min_eqs_for_a_cluster = 30; mult_incr = 1;
% Syntax: [n0,analy, value_counts] = clustering_analysis...
%                 (xs, ys, zs, strike_incr, dip_incr, width, mult_incr, min_eqs_for_a_cluster)
% 
%
%The solution is not sensitive to angles greater than 180. 181 deg is the same as 1 deg
nhypos = length(xs);
strikes = 0:strike_incr:175; %pos is anticlockwise. 
dips = 0:dip_incr:90;

% Initializing array: analy = [1 2 3 - x, y, z
%                              4,5,6 - rxp, ryp, rzp
%                              7,8   - neqs_in_window, st (dist at max eqs)
%                              9, 10, 11] - strike, dip, unique_value

analy = [xs' ys' zs' zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1) ...
    zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1) zeros(nhypos,1)  zeros(nhypos,1)]; 

unique_value = 0;

% Analyzing for each combination of strike and dip. We do not need the rake.
for j = 1: length(strikes)
    for jj = 1:length(dips)
        
        % Getting a unique value
        unique_value = unique_value + 1; 
        
        con=pi/180.;
        strike=strikes(j).*con;
        dip=dips(jj).*con;

        % extracting the xs, ys and zs because their relative location will change
        % during the analysis.
        xs_now = analy(:,1)';
        ys_now = analy(:,2)';
        zs_now = analy(:,3)';

        % removing the mean for the rotation to be about origins. 
        xs_now = xs_now - mean(xs_now); 
        ys_now = ys_now - mean(ys_now);
        zs_now = zs_now - mean(zs_now);

        % concatenating x,y and z
        R=[xs_now ; ys_now ;zs_now];

        %************** rotate into strike direction **********************
        %Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];
        Dstrike=[ cos(strike)  -sin(strike) 0 ; sin(strike) cos(strike) 0 ; 0 0 1];

        Rstrike=Dstrike*R;

        rxp(1:nhypos) = Rstrike(1,1:nhypos);
        ryp(1:nhypos) = Rstrike(2,1:nhypos);
        rzp(1:nhypos) = Rstrike(3,1:nhypos);

        %************** rotate into dip direction ********************
        Rstrike(1,1:nhypos) = Rstrike(1,1:nhypos) - mean(Rstrike(1,1:nhypos));
        Rstrike(2,1:nhypos) = Rstrike(2,1:nhypos) - mean(Rstrike(2,1:nhypos));
        Rstrike(3,1:nhypos) = Rstrike(3,1:nhypos) - mean(Rstrike(3,1:nhypos));

        % Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];
        Ddip=[ cos(dip) 0 -sin(dip); 0 1  0 ; sin(dip) 0 cos(dip)];

        Rdip=Ddip*Rstrike;

        % add bac the barycenter
        Rdip(1,1:nhypos) = Rdip(1,1:nhypos) + mean(Rstrike(1,1:nhypos));
        Rdip(2,1:nhypos) = Rdip(2,1:nhypos) + mean(Rstrike(2,1:nhypos));
        Rdip(3,1:nhypos) = Rdip(3,1:nhypos) + mean(Rstrike(3,1:nhypos));

        rxp(1:nhypos) = Rdip(1,1:nhypos);
        ryp(1:nhypos) = Rdip(2,1:nhypos);
        rzp(1:nhypos) = Rdip(3,1:nhypos);

        % Take the depth to have a minimum of zero. Thi does not overright the
        % original dataset.
        rzp = rzp-min(rzp);
        
        % Concatenate these result with the original dataet to kep track of each
        % point. Sort the data based on the new depths.
        analy(:,4) = rxp';
        analy(:,5) = ryp';
        analy(:,6) = rzp';

        % Sorting in the new depth
        analy = sortrows(analy,6);

        % Stepwisely move through the data in depth and the number of EQs in each
        % block. we need to specify the width of the block. 
        st_array = 0:max(rzp+1); 

        for i = 1:length(st_array)
            st = st_array(i);
            
            % Determine the no of eqs in each depth window
            index = find(analy(:,6)>=st & analy(:,6)<(st+width));
            neqs_in_window = length(index);

            % For each eq in a data window, we want the no of eqs in the
            % cluster, dist index, strike, dip and a unique value against 
            % the hypocenter but before that, we wnt to
            % check if the no of eqs assigned to thi hypocenter is less
            % than the current value. Also, if the no of eqss is >=
            % min_eqs_for_a_cluster specified by the user.
            for ii = 1: neqs_in_window
                if (analy(index(ii),7) < neqs_in_window) && ...
                        (neqs_in_window >= min_eqs_for_a_cluster)
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

% When the analysis of strike and dips are completed, sort analy by the
% unique value. Determine the unique values in the unique_value and count
% the no of eqs for each unique value. Sort in descending no of eqs.
% Clusters with the same strike and dip will be differentiated later by the
% dist (i.e. st) in the depth window.
analy = sortrows(analy,11);

[clus,ia,ic] = unique(analy(:,11));
a_counts = accumarray(ic,1);
value_counts = [clus, a_counts];
value_counts = sortrows(value_counts,2,'descend');

% For another matrix value_counts = [unique_value no_of_eqs dist_index
% strike dip]
for i= 1:length(value_counts(:,1))
    dist_index = analy(analy(:,11)== value_counts(i,1),8); value_counts(i,3) = dist_index(1);
    strike_a = analy(analy(:,11)== value_counts(i,1),9); value_counts(i,4) = strike_a(1);
    dip_a = analy(analy(:,11)== value_counts(i,1),10); value_counts(i,5) = dip_a(1);

end

value_counts = round(value_counts);

value_counts(:,6) = 0;

% Combining clusters with similar strike, dip, dist_index
% dstrike <= mult_incr * strike_incr
% ddip <= mult_incr * dip_incr
% ddist_index <= 1
% value_counts(:,6) == 0 i.e. if it doesnot have any previously assigned
% value to avoid overwritten.
for i= 1:length(value_counts(:,1))
    
    if value_counts(i,6) == 0
        
        value_counts((abs(value_counts(:,3) - value_counts(i,3))<=1) & ...
            (abs(value_counts(:,4) - value_counts(i,4))<=mult_incr*strike_incr) & ...
            (abs(value_counts(:,5) - value_counts(i,5))<=mult_incr*dip_incr)& ...
            (value_counts(:,6) == 0) ,6) = i; % 

    end
end

% sorting again and based on the unique value from value_counts(:,6),
% overwrite the strike, dip, dist_index and the sum of the eqs in the analy
% matrix
uni = unique(value_counts(:,6));

for i=1:length(uni)
    uniq_index = value_counts(value_counts(:,6) == uni(i),1);

    for j = 1:length(uniq_index)
        
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

% Determine the no of clusters with eqs >= min_eqs_for_a_cluster
n0 = length(value_counts(value_counts(:,2)>=min_eqs_for_a_cluster));

% ------------------------------------------------------------------------
% Creating figures
% ------------------------------------------------------------------------
fig = figure;
ax1 = subplot(1,2,1); %ax1 = 
plot3(xs,ys,zs,'o');
axis equal; title('Input Hypocenters');
xlabel('X km'); ylabel('Y km'); zlabel('Z km'); grid on

ax2 = subplot(1,2,2); % ax2 = 
for ncluster = 1:length(cluster)
    
    if value_counts(ncluster,2) >= min_eqs_for_a_cluster
        xsc = analy(analy(:,11) == cluster(ncluster),1);
        ysc = analy(analy(:,11) == cluster(ncluster),2);
        zsc = analy(analy(:,11) == cluster(ncluster),3);

        plot3(xsc,ysc,zsc,'o'); hold on;
        
    else
        xsc = analy(analy(:,11) == cluster(ncluster),1);
        ysc = analy(analy(:,11) == cluster(ncluster),2);
        zsc = analy(analy(:,11) == cluster(ncluster),3);

        plot3(xsc,ysc,zsc,'ko'); hold on;
    
    end
end

axis equal; title('Clustered Data'); 
xlabel('X km'); ylabel('Y km'); zlabel('Z km'); grid on

analy = [analy(:,1) analy(:,2) analy(:,3) analy(:,7) analy(:,8) analy(:,9) analy(:,10) analy(:,11)]; 


hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); 
rotate3d on



% OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)

%toc 

%*************************** END ******************************************