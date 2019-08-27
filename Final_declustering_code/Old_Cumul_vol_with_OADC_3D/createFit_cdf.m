function pd1 = createFit_cdf(a)
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

% This function was automatically generated on 14-May-2019 17:11:53

% Output fitted probablility distribution: PD1

% Data from dataset "a data":
%    Y = a

% Force all inputs to be column vectors
a = a(:);

% Prepare figure
%clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "a data"
[CdfY,CdfX] = ecdf(a,'Function','cdf');  % compute empirical function
hLine = stairs(CdfX,CdfY,'Color',[0.333333 0 0.666667],'LineStyle','-', 'LineWidth',1);
xlabel('Movement (km)');
ylabel('Cumulative probability')
LegHandles(end+1) = hLine;
LegText{end+1} = 'data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('gamma',[ 3.042390831227, 0.1977322266891])
pd1 = fitdist(a, 'normal');%'gamma' 'gamma'
YPlot = cdf(pd1,XGrid);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'fit 1';

% Adjust figure
box on;
grid on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
set(hLegend,'Interpreter','none');
