function plot_figures_for_OADC(nn)

%load('Final_result_Simul.1.OADC_results/Simul.1.OADC.saved_variables.mat')

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global Strike Dip FM_file dist2FM_threshold
global xb_tmp_i yb_tmp_i zb_tmp_i
global xv_tmp_i yv_tmp_i zv_tmp_i
global xt_tmp_i yt_tmp_i zt_tmp_i
global vec_plane_tmp_i
global Nt_tmp_i lambda3_tmp_i
global L_tmp_i W_tmp_i Strike_tmp_i Dip_tmp_i
global index use_glo_var con_tol Kfaults

n=nn;

% load up arrays with good cluster parameters    
xb = xb_tmp_i(:,:,index(n)); yb = yb_tmp_i(:,:,index(n)); 
zb=zb_tmp_i(:,:,index(n)); % Barycenters
xv = xv_tmp_i(:,:,index(n)); yv=yv_tmp_i(:,:,index(n)); 
zv=zv_tmp_i(:,:,index(n)); % fault plane vertices
xt=xt_tmp_i(:,:,index(n)); yt=yt_tmp_i(:,:,index(n)); 
zt=zt_tmp_i(:,:,index(n)); % hypocenter location in a cluster
vec_plane=vec_plane_tmp_i(:,:,index(n)); % eigenvector that describes each plane
Nt=Nt_tmp_i(:,:,index(n)); % number of events in each trial cluster
lambda3=lambda3_tmp_i(:,:,index(n)); % minimum eigenvalue
L=L_tmp_i(:,:,index(n)); W=W_tmp_i(:,:,index(n)); 
Strike=Strike_tmp_i(:,:,index(n)); Dip=Dip_tmp_i(:,:,index(n)); % fault plane parameters
            
if use_glo_var == 1

    JFINAL=faultcluster(con_tol,Kfaults);
else

    JFINAL=Copy_of_faultcluster(con_tol,Kfaults);
end

Strike
Dip

lambda3
JFINAL

figure;                
picname=['Final Model ' int2str(n)]; 
simul_tag= 'Simul.1.OADC'; 
%n = 2; xv_disp = xv_tmp_i(:,:,index(n));
%yv_disp = yv_tmp_i(:,:,index(n)); zv_disp = zv_tmp_i(:,:,index(n));
datplot_with_colors(xs,ys,zs,Kfaults,xv,yv,zv,picname,simul_tag);
