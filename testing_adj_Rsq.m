clear all; close all; clc;

%Create two variables, x and y, from the first two columns of the count 
% variable in the data file count.dat:


x = 0:10;
y = 5*x + 3; 

% Perturbing y
y(5) = y(5) + 2;
y(4) = y(4) + 8;

% Call polyfit to generate a cubic fit to predict y from x:

p = polyfit(x,y,1);

%Call polyval to use the coefficients in p to predict y, naming the result yfit:

yfit = polyval(p,x);

%Compute the residual values as a vector of signed numbers:

yresid = y - yfit;

% Square the residuals and total them to obtain the residual sum of squares:

SSresid = sum(yresid.^2);

%Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:

SStotal = (length(y)-1) * var(y);

%Compute simple R2 for the cubic fit using the formula given in the introduction of this topic:

rsq = 1 - SSresid/SStotal

%Finally, compute adjusted R2 to account for degrees of freedom:

rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p))
