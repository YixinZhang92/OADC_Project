function pd1 = createFit_pf(a)
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(A)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 14-May-2019 17:36:59

% Output fitted probablility distribution: PD1

% Data from dataset "a data":
%    Y = a

% Force all inputs to be column vectors
a = a(:);

% Prepare figure
%clf;
hold on;
LegHandles = []; LegText = {};
CustomDist = internal.stats.dfgetdistributions('gamma');
probplot({CustomDist,[ 3.042390831227, 0.1977322266891]});
title('');


% --- Plot data originally in dataset "a data"
hLine = probplot(gca,a,[],[],'noref'); % add data to existing plot
set(hLine,'Color',[0.333333 0 0.666667],'Marker','o', 'MarkerSize',6);
xlabel('Movement (km)');
ylabel('Probability')
LegHandles(end+1) = hLine;
LegText{end+1} = 'a data';


% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('gamma',[ 3.042390831227, 0.1977322266891])
pd1 = fitdist(a, 'normal');%'gamma'
hLine = probplot(gca,pd1);
set(hLine,'Color',[1 0 0],'LineStyle','-', 'LineWidth',2);
LegHandles(end+1) = hLine;
LegText{end+1} = 'fit 1';

% Adjust figure
box on;
grid on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
set(hLegend,'Interpreter','none');
