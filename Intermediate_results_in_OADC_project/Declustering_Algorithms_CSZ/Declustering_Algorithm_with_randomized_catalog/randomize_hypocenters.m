function [x_rand, y_rand, z_rand]=randomize_hypocenters(x,y,z,rand_rang)

N = length(x);

random_x = (rand(N,1)*2*rand_rang) - rand_rang;
random_y = (rand(N,1)*2*rand_rang) - rand_rang;

x_rand = x + random_x;
y_rand = y + random_y;
z_rand = z;