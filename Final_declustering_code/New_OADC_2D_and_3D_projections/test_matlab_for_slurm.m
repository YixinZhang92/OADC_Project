function test_matlab_for_slurm()
%clear all; close all; clc;

tic 
simul_tag  = 'test';

x = 1:10;
y = 1:10;
z = 1:10;

% Ploting figures
Fig1 = figure('Name','Clustered hypocenters with fault model','Position', get(0, 'Screensize'));

plot3(x,y,z,'ro')

grid on; axis equal;
title('Fault Model (Constraint: Lambda2 only)');
set(gca, 'fontsize', 18); shg

% saving all variables to file
savevar_filename = [simul_tag '.saved_variables.mat'];
save(savevar_filename)

% Printing figure to file
fig_filename = [simul_tag '.faultmodel.lambda2_only.png'];
F1    = getframe(Fig1);
imwrite(F1.cdata, fig_filename, 'png')
savefig(Fig1,[fig_filename(1:end-4) '.fig']);


