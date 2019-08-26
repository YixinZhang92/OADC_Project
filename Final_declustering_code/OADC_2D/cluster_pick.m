function varargout = cluster_pick(varargin)
% CLUSTER_PICK MATLAB code for cluster_pick.fig
%      CLUSTER_PICK, by itself, creates a new CLUSTER_PICK or raises the existing
%      singleton*.
%
%      H = CLUSTER_PICK returns the handle to a new CLUSTER_PICK or the handle to
%      the existing singleton*.
%
%      CLUSTER_PICK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_PICK.M with the given input arguments.
%
%      CLUSTER_PICK('Property','Value',...) creates a new CLUSTER_PICK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_pick_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_pick_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_pick

% Last Modified by GUIDE v2.5 02-Oct-2015 10:11:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_pick_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_pick_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cluster_pick is made visible.
function cluster_pick_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_pick (see VARARGIN)

% Choose default command line output for cluster_pick
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_pick wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_pick_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function azimuth_slider_Callback(hObject, eventdata, handles)
% hObject    handle to azimuth_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Determine the azimuth increment/decrement to use in the rotation

azimuth=get(handles.azimuth_slider,'Value');
%fprintf('azimuth= %g\n',azimuth);

set(handles.azimuth_edit,'String',num2str(azimuth));



% --- Executes during object creation, after setting all properties.
function azimuth_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to azimuth_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function elevation_slider_Callback(hObject, eventdata, handles)
% hObject    handle to elevation_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

elevation=get(handles.elevation_slider,'Value');

set(handles.elevation_edit,'String',num2str(elevation));


% --- Executes during object creation, after setting all properties.
function elevation_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elevation_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function azimuth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to azimuth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of azimuth_edit as text
%        str2double(get(hObject,'String')) returns contents of azimuth_edit as a double


% --- Executes during object creation, after setting all properties.
function azimuth_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to azimuth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elevation_edit_Callback(hObject, eventdata, handles)
% hObject    handle to elevation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elevation_edit as text
%        str2double(get(hObject,'String')) returns contents of elevation_edit as a double


% --- Executes during object creation, after setting all properties.
function elevation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elevation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pick_pushbutton.
function pick_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pick_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Pick a polygon to delimit a fault seismicity cluster, subset the
%  seismicity, calculate a rectangular fault plane model using principal
%  components analysis

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

h1=handles.density_plot_axes;

%  Pick the polygon
[x_poly,y_poly]=getline(h1,'closed');

%  subset the seismicity within the polygon
% determine if a hypocenter is inside or outside the polygon
IN=inpolygon(xd,yd,x_poly,y_poly);

% subset the seismicity in this polygon
% each transformed event is in the same order as the input data
k=0;
for kk=1:N;
    
    if IN(kk) == 1;
        k=k+1;
        xn(k)=xs(kk);
        yn(k)=ys(kk);
        zn(k)=zs(kk);
    end;
    
end;

%  save the clustered seismicity
if k >= 0;
    Cluster_number=Cluster_number+1;
    N_cluster(Cluster_number)=k;
    xclus(Cluster_number,1:N_cluster(Cluster_number))=xn(1:N_cluster(Cluster_number));
    yclus(Cluster_number,1:N_cluster(Cluster_number))=yn(1:N_cluster(Cluster_number));
    zclus(Cluster_number,1:N_cluster(Cluster_number))=zn(1:N_cluster(Cluster_number));  
end

%  check plot
%figure;
%plot3(xclus(Cluster_number,:),yclus(Cluster_number,:),zclus(Cluster_number,:),'ro');

%**************************************************************************
%  Principal Components analysis
%  first compute the Barycenter of this cluster
kk=Cluster_number;
nclus=N_cluster(Cluster_number);
xb(kk)=mean(xclus(kk,1:nclus));
yb(kk)=mean(yclus(kk,1:nclus));
zb(kk)=mean(zclus(kk,1:nclus));

% compute the covariance matrix for this cluster
Cxy=cov( [xclus(kk,1:nclus)' yclus(kk,1:nclus)' zclus(kk,1:nclus)'],0);
        
% compute the eigenvalues and eigenvectors for this cluster
[V,De]=eig(Cxy);
            
% calculate fault plane parameters from the eigen results
% and calculate the vertices of the fault plane
[L(kk),W(kk),Strike(kk),Dip(kk),xv(kk,:),yv(kk,:),zv(kk,:)] = fltplane(V,De,xb(kk),yb(kk),zb(kk));

% save the plane unit normal vector and eigenvalue
vec_plane(kk,1:3)=V(1:3,1);
lambda3(kk)=sqrt(12.*De(1,1));

% another checkplot
picname='cluster check';
datplot2(xn,yn,zn,kk,xv,yv,zv,picname);


%  Now plot the density plot on the GUI plot axes
PLOT_Latest=1;
h1=handles.density_plot_axes;
cla(h1,'reset');
set(h1,'NextPlot','add');

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
elevation=str2num(get(handles.elevation_edit,'String'));
azimuth=str2num(get(handles.azimuth_edit,'String'));
zfactor=get(handles.zoom_slider,'Value');

view_flag=0;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');



% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plot_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plot_density_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_density_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot the current density plot window in a separate figure

figure('Name','Density Plot');
h3=axes;

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'))
elevation=str2num(get(handles.elevation_edit,'String'))
azimuth=str2num(get(handles.azimuth_edit,'String'))
zfactor=get(handles.zoom_slider,'Value')

view_flag=0

plot_density(h3,psi_angle,elevation,azimuth,zfactor,view_flag);




% --------------------------------------------------------------------
function plot_faults_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_faults_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot all seismicity and faults in a separate figure

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

% plot a rendition of the data and the planes that fit the data
figure('Name','All seismicity and Faults');

hold on;

plot3(xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);

for k=1:Cluster_number;
    fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end;

hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
grid on;
view(3);

% --------------------------------------------------------------------
function read_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to read_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

%  Initialize all working space
initialize_all;

%  Read a seismicity file consisting of (x,y,z), depth is negative

%   Call uigetfile menu

[filename,dirname]=uigetfile('*.*');
if filename == 0; return;end;

infile=strcat(dirname,filename);
fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data1=data;
data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

% plot the input data on a separate figure
% 
figure('Name','Seismicity');
h0=axes;
plot3(h0,xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
axis equal;
grid on;
title(h0,'Input hypocenters');
xlabel(h0,'X km');
ylabel(h0,'Y km');
zlabel(h0,'Z km');

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
elevation=get(handles.elevation_slider,'Value');
azimuth=get(handles.azimuth_slider,'Value');
zfactor=get(handles.zoom_slider,'Value');

view_flag=0;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');

% --------------------------------------------------------------------
function write_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to write_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Write the clustered seismicity and fault plane models to a file

global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

%  Format of clustered seismicity
%       Cluster_number - index number of cluster
%       N_cluster - number of events in cluster
%       xclus(Cluster_number,1:N_cluster(Cluster_number)) - event x coordinate
%       yclus(Cluster_number,1:N_cluster(Cluster_number)) - event y coordinate
%       zclus(Cluster_number,1:N_cluster(Cluster_number)) - event z coordinate

%  Fault Plane Model (k=Cluster_number)
%       L(k)  - length
%       W(k)  - width
%       Strike(k)
%       Dip(k)
%       xv(k,1:4)    Rectangular fault plane vertices
%       yv(k,1:4)
%       zv(k,1:4)
%       lambda3(k)   fault plane eigenvalue

%  uiputfile dialog

[filename,dirname]=uiputfile('*.*');
if filename == 0; return;end;
    
faultfile=strcat(dirname,filename);

fid=fopen(faultfile,'w');

fprintf(fid,'%s\n',[faultfile]);

for k=1:Cluster_number;
    
    %  cluster index
    fprintf(fid,'%i \n',k);
    
    %  fault geometry
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g %12.5g \n',[L(k) W(k) Strike(k) Dip(k) lambda3(k)]);
    
    %  rectangular fault vertices
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g \n',[xv(k,1:4)]);
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g \n',[yv(k,1:4)]);
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g \n',[zv(k,1:4)]);
    
    %  cluster seismicity
    n=N_cluster(k);
    fprintf(fid,'%i \n',n);
    
    for kk=1:n
        fprintf(fid,'%12.5g %12.5g %12.5g \n',[xclus(k,kk) yclus(k,kk) zclus(k,kk)]);
    end
    
end;

fclose(fid);
        


function zoom_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_edit as text
%        str2double(get(hObject,'String')) returns contents of zoom_edit as a double


% --- Executes during object creation, after setting all properties.
function zoom_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function zoom_slider_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
PLOT_Latest=0;

zfactor=get(handles.zoom_slider,'Value');
%fprintf('elevation= %g\n',elevation);

set(handles.zoom_edit,'String',num2str(zfactor));

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
elevation=get(handles.elevation_slider,'Value');
azimuth=get(handles.azimuth_slider,'Value');

view_flag=0;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');

% --- Executes during object creation, after setting all properties.
function zoom_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in translate_pushbutton.
function translate_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to translate_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D

%  Pick a new (x,y) relative origin for data

h1=handles.density_plot_axes;
set(h1,'Nextplot','add');

%  Pick a new origin point
[x_test,y_test]=ginput(1);

%  convert into original coordinates using the inverse rotation matrix
xm_old=[x_test ; y_test ; 0.0]
xm_new=DI*xm_old;

x_orig=xm_new(1) + x_orig;
y_orig=xm_new(2) + y_orig;

n_orig=n_orig+1;

%  Now plot the density plot on the GUI plot axes
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
elevation=get(handles.elevation_slider,'Value');
azimuth=get(handles.azimuth_slider,'Value');
zfactor=get(handles.zoom_slider,'Value');

view_flag=0;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');

return


% --- Executes on button press in refresh_pushbutton.
function refresh_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  refresh seismicity plot
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D

D=[ 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
x_orig=0.0;
y_orig=0.0;

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
elevation=get(handles.elevation_slider,'Value');
azimuth=get(handles.azimuth_slider,'Value');
zfactor=get(handles.zoom_slider,'Value');

view_flag=0;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');

% --------------------------------------------------------------------
function plot_faults_only_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_faults_only_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot all faults in a separate figure

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

figure('Name','Faults Only');

hold on;

for k=1:Cluster_number;
    fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end;

hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
grid on;
view(3);

% --------------------------------------------------------------------
function plot_clusters_faults_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_clusters_faults_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot cluster seismicity with faults in a separate figure

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

% plot a rendition of the data and the planes that fit the data
figure('Name','Cluster seismicity and Faults');

hold on;

for k=1:Cluster_number;
   
    kc=N_cluster(k);
    plot3(xclus(k,1:kc),yclus(k,1:kc),zclus(k,1:kc),'o','MarkerEdgeColor','k','MarkerFaceColor','k');
    fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end;

hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
grid on;
view(3);

% --- Executes on button press in apply_azimuth_pushbutton.
function apply_azimuth_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to apply_azimuth_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global elevation_old azimuth_old

PLOT_Latest=0;

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
azimuth=str2num(get(handles.azimuth_edit,'String'));
elevation=str2num(get(handles.elevation_edit,'String'));

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

zfactor=get(handles.zoom_slider,'Value');

view_flag=3;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');



% --- Executes on button press in apply_elevation_pushbutton.
function apply_elevation_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to apply_elevation_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global elevation_old azimuth_old
PLOT_Latest=0;

elevation=str2num(get(handles.elevation_edit,'String'));

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));
azimuth=str2num(get(handles.azimuth_slider,'String'));
zfactor=get(handles.zoom_slider,'Value');

view_flag=2;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');


% --- Executes on slider movement.
function rotate_about_x_slider_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_about_x_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

psi_angle=get(handles.rotate_about_x_slider,'Value');

set(handles.rotate_about_x_edit,'String',num2str(psi_angle));

% --- Executes during object creation, after setting all properties.
function rotate_about_x_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_about_x_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function rotate_about_x_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_about_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotate_about_x_edit as text
%        str2double(get(hObject,'String')) returns contents of rotate_about_x_edit as a double


% --- Executes during object creation, after setting all properties.
function rotate_about_x_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_about_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply_x_rotation_pushbutton.
function apply_x_rotation_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to apply_x_rotation_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global elevation_old azimuth_old

%  rotate about the x axis

PLOT_Latest=0;

azimuth=str2num(get(handles.azimuth_edit,'String'));
elevation=str2num(get(handles.elevation_edit,'String'));
psi_angle=str2num(get(handles.rotate_about_x_edit,'String'));

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

zfactor=get(handles.zoom_slider,'Value');

view_flag=1;

plot_density(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');


% --- Executes on key press with focus on apply_x_rotation_pushbutton and none of its controls.
function apply_x_rotation_pushbutton_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to apply_x_rotation_pushbutton (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
