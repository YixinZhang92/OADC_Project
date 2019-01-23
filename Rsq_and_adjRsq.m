function [Rsq, adj_Rsq] = Rsq_and_adjRsq(x, y, d, st)
% clear all; close all; clc;
%Create two variables, x and y, from the first two columns of the count 
%  x= 0:10;
% y = 5*x + 4;
% d =1;

% Call polyfit to generate a cubic fit to predict y from x:
p = polyfit(x,y,d);

if nargin < 4 % i.e. when the constant yfit is not given
    %Call polyval to use the coefficients in p to predict y, naming the result yfit:
    yfit = polyval(p,x);
    
else
    yfit = st*ones(length(x),1)';
end

%Compute the residual values as a vector of signed numbers:
yresid = y - yfit;

% Square the residuals and total them to obtain the residual sum of squares:
SSresid = sum(yresid.^2);

%Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(y)-1) * var(y);

%Compute simple R2 for the cubic fit using the formula given in the introduction of this topic:
Rsq = 1 - SSresid/SStotal;

%Finally, compute adjusted R2 to account for degrees of freedom:
adj_Rsq = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p));


% function c = addme(a,b)
%     switch nargin
%         case 2
%             c = a + b;
%         case 1
%             c = a + a;
%         otherwise
%             c = 0;
%     end
% end