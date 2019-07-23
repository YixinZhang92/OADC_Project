function mat = vect2mat_manual(vec,matcol)
% clear all; close all; clc;
% vec= 1:12;
% matcol = 4;

N=length(vec);
vec = reshape(vec,[1,N]);

no_rows = ceil(N/matcol);
no_el_needed = no_rows*matcol;

n_padded = no_el_needed- N;
zero_pad = zeros(1,n_padded);

vec_now = [vec zero_pad];
mat = reshape(vec_now,[matcol no_rows])';