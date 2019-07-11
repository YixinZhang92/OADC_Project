function varargout = focal(varargin)
% FOCAL M-file for focal.fig
%      FOCAL, by itself, creates a new FOCAL or raises the existing
%      singleton*.
%
%      H = FOCAL returns the handle to a new FOCAL or the handle to
%      the existing singleton*.
%
%      FOCAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOCAL.M with the given input arguments.
%
%      FOCAL('Property','Value',...) creates a new FOCAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before focal_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to focal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help focal

% Last Modified by GUIDE v2.5 20-Sep-2003 14:57:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @focal_OpeningFcn, ...
                   'gui_OutputFcn',  @focal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before focal is made visible.
function focal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to focal (see VARARGIN)

% Choose default command line output for focal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%set(handles.focal_tool,'MenuBar','figure');

global DATA_READ_FLAG

DATA_READ_FLAG=0;       %  Initialize data file reading flag
%
% -------------------------------------------------------------------------
% UIWAIT makes focal wait for user response (see UIRESUME)
% uiwait(handles.focal_tool);

% --- Executes on button press in read_arrival_data.
function read_arrival_data_Callback(hObject, eventdata, handles)
% hObject    handle to read_arrival_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function creates polarity and amplitude data structure from 
% either synthetic or real data.
% 
% Outputs: p_pol, sv_pol, sh_pol [azimuth, Amplitude and incidence angle] 
%          vp, vs, dens at the source depth
%
% -------------------------------------------------------------------------
%  set global variable for arrivals file
global p_pol sv_pol sh_pol vp vs dens %strike_incr dip_incr rake_incr  
% This is because the initial incr may be changed.
global dist az depth halfspace vel_model path_to_SAC path_to_travelt
global w_pol w_ratio use_sv_pol use_sv_amp depth_error z_displ r_displ t_displ
%
% Run the parameters_in.m file to get some of the input parameters for focal.
parameters_in()
%
% -------------------------------------------------------------------------
if synthetic == 1         
    % Determining synthetic polarity data using the given info. 
    % The weights of all amplitudes are automatically 5.
    [p_pol, sv_pol, sh_pol, vp, vs, dens] = synth_polarities...
        (dist,az,depth, halfspace, vel_model, strike, dip, rake,filename,...
         path_to_SAC, path_to_travelt, z_displ, r_displ, t_displ);
  
else
    % Real dataset  
    [p_pol, sv_pol, sh_pol, vp, vs, dens] = real_polarities...
        (dist,az,depth, halfspace, vel_model, p_data, sv_data, sh_data,...
         path_to_SAC, path_to_travelt);
 
end

p_pol
sv_pol
sh_pol
%
% -------------------------------------------------------------------------
%   send message to message strip
if synthetic == 1  
    set(handles.message_strip,'String',['Reading in file ',filename,...
        ' and setting global variables']);
else
    set(handles.message_strip,'String',...
        'Reading in real data and setting global variables');
end
%
% -------------------------------------------------------------------------



% --- Executes on button press in perform_grid_search.
function perform_grid_search_Callback(hObject, eventdata, handles)
% hObject    handle to perform_grid_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%   This function calculates the following:
%         1. The P, SV, and SH wave nodal surfaces that fits the input wave 
%            polarity and |SV|/SH, SH/P and |SV|/P amplitude ratios data
%            using grid search method,
%         3. Plot P-wave nodal surfaces of the results if the user says YES
%            in PLOT NODAL SURFACES. and the Swave polarization vector
%            using the SV and SH amplitude data.
%
% -------------------------------------------------------------------------
%  extracting global variables
global p_pol sv_pol sh_pol vp vs dens strike_incr dip_incr rake_incr 
global grid_search_results w_pol w_ratio use_sv_amp use_sv_pol 
global z_displ r_displ t_displ
%
% -------------------------------------------------------------------------
%   send message to message strip
set(handles.message_strip,'String','Grid Search in progress ...');
%
% -------------------------------------------------------------------------
strike_incr = str2double(get(handles.strike_inc,'String'));
dip_incr = str2double(get(handles.dip_inc,'String'));
rake_incr = str2double(get(handles.rake_inc,'String'));
% -------------------------------------------------------------------------
% Performing grid search method to determine best-fit focal mechanisms
grid_search_results = vectorized_grid_search_analyses...
      (strike_incr,dip_incr,rake_incr,p_pol,sv_pol,sh_pol,vp,vs,dens, ...
       w_pol, w_ratio, use_sv_amp, use_sv_pol, z_displ, r_displ, t_displ);
% 
% -------------------------------------------------------------------------
%   send message to message strip 
strike1 = grid_search_results(1,1);
dip1 = grid_search_results(1,2);
rake1 = grid_search_results(1,3);

set(handles.message_strip,'String',strcat...
    ('Correct FM: Strike: ',num2str(strike1),' ; Dip: ',...
    num2str(dip1),' ; Rake: ',num2str(rake1)));
%
% -------------------------------------------------------------------------
% send message to number of focal mechanisms found.
[numr, ~] = size(grid_search_results);
set(handles.number_mech,'String',num2str(numr));
%
% -------------------------------------------------------------------------
% Plots
% ------------------------------------------------------------------------- 
% Plotting the PTI axes of the two best sets of FMs on focalsphere
axes(handles.pwavenodes); cla; set(handles.pwavenodes,'NextPlot','add',...
    'XTick',[],'YTick',[],'ZTick',[]);
    plot_focal_figures(grid_search_results, 'p_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens, z_displ, r_displ, t_displ)
set(handles.pwavenodes,'NextPlot','replace');
%
% -------------------------------------------------------------------------  
% Plotting the P polarities & P-nodal planes of the best FMs on focalsphere
axes(handles.svwavenodes); cla;set(handles.svwavenodes,'NextPlot','add',...
    'XTick',[],'YTick',[],'ZTick',[]);    
    plot_focal_figures(grid_search_results, 'sv_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens, z_displ, r_displ, t_displ)
set(handles.svwavenodes,'NextPlot','replace'); 
%
% ------------------------------------------------------------------------- 
% Plotting the Swave polarization vector and P-nodal planes of the 
% best FMs on focalsphere
axes(handles.shwavenodes); cla;set(handles.shwavenodes,'NextPlot','add',...
    'XTick',[],'YTick',[],'ZTick',[]);    
    plot_focal_figures(grid_search_results, 'sh_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens, z_displ, r_displ, t_displ)
set(handles.shwavenodes,'NextPlot','replace');
%
% -------------------------------------------------------------------------
fig10 = figure (10);
clf 
   plot_focal_figures(grid_search_results, 'p_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens, z_displ, r_displ, t_displ)
   axis off;
   print(fig10, 'fig_PTI_axes.png', '-dpng') 
% -------------------------------------------------------------------------
fig11 = figure (11);
clf 
   plot_focal_figures(grid_search_results, 'sv_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens, z_displ, r_displ, t_displ)
   axis off;
   print(fig11, 'fig_P_pol.png', '-dpng') 
% -------------------------------------------------------------------------
fig12 = figure (12);
clf 
   plot_focal_figures(grid_search_results, 'sh_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens, z_displ, r_displ, t_displ)
   axis off;
   print(fig12, 'fig_Swave_polariz.png', '-dpng') 
% -------------------------------------------------------------------------  
%save grid_search_results.mat grid_search_results
% [strike2,dip2,rake2] = calc_aux_plane(strike1,dip1,rake1);
% 
% %   send message to message strip 
% set(handles.message_strip,'String',strcat...
%     ('Correct FM: Strike: ',num2str(strike1),' ; Dip: ',...
%     num2str(dip1),' ; Rake: ', num2str(rake1),' (Aux Plane= Strike: ',...
%     num2str(strike2),' ; Dip:  ',num2str(dip2),' ; Rake:  ',num2str(rake2),')'));
% -------------------------------------------------------------------------



% --- Executes on button press in write_mechanisms.
function write_mechanisms_Callback(hObject, eventdata, handles)
% hObject    handle to write_mechanisms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
global grid_search_results
%
filename = 'write_mechanisms.txt';
fid = fopen(filename,'w'); % open the original vtk file 
%
% test if the file exists and able to open
if(fid==-1) 
    error('Can''t open the file.');
end
%
dlmwrite(filename,grid_search_results,'delimiter','\t')
%
%   send message to message strip
set(handles.message_strip,'String',...
    'Writing only the best Mechanisms to file');
              
         
% --- Executes on button press in read_mechanisms.
function read_mechanisms_Callback(hObject, eventdata, handles)
% hObject    handle to read_mechanisms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% -------------------------------------------------------------------------            
global p_pol sv_pol sh_pol strike_incr dip_incr rake_incr use_sv_pol
global w_pol w_ratio dist az depth halfspace vel_model depth_error use_sv_amp
global path_to_SAC path_to_travelt z_displ r_displ t_displ
%
% -------------------------------------------------------------------------
%   send message to message strip
set(handles.message_strip,'String',...
    'Performing error analyses using errors in polarities and focal depths');
%
% -------------------------------------------------------------------------
% Perform error analysis by changing the polarities at each station
% focal_error_analysis(strike_incr,dip_incr,rake_incr,p_pol,sv_pol,...
% sh_pol,vp,vs,dens,error_analysis)
% focal_error_analysis(strike_incr,dip_incr,rake_incr,...
%    p_pol,sv_pol,sh_pol,vp,vs,dens, w_pol, w_ratio, use_sv_amp, use_sv_pol)

% perform error analysis by changing the focal depths and change the
% polarities at each station at each iteration of depth.
focal_error_analysis_focaldepth(strike_incr,dip_incr,rake_incr,...
    p_pol,sv_pol,sh_pol, w_pol, w_ratio, dist, az, depth, halfspace, ...
    vel_model, depth_error, use_sv_amp, use_sv_pol, path_to_SAC, ...
    path_to_travelt, z_displ, r_displ, t_displ)
% -------------------------------------------------------------------------


%**************************************************************************
% --- Executes on button press in calc_nodal_surf.
function calc_nodal_surf_Callback(hObject, eventdata, handles)
% hObject    handle to calc_nodal_surf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%   This function calculates the P, SV, and SH wave nodal surfaces
% -------------------------------------------------------------------------
%  set global variable for arrivals file
global vp vs dens nodeflag axflag dataflag stype mom
global p_pol sv_pol sh_pol nmech
global z_displ r_displ t_displ
%  
% default value of nmech
nmech=1;
%
% -------------------------------------------------------------------------
% % converting P-wave amplitudes to polarities
% p_pol(:,1) = p_pol(:,1)./abs(p_pol(:,1));
% sv_pol(:,1) = sv_pol(:,1)./abs(sv_pol(:,1));
% sh_pol(:,1) = sh_pol(:,1)./abs(sh_pol(:,1));
%
% -------------------------------------------------------------------------
%   Get source velocity and density structure
vp=str2num(get(handles.vp_input,'String'));
vs=str2num(get(handles.vs_input,'String'));
dens=str2num(get(handles.dens_input,'String'));
%
% -------------------------------------------------------------------------
% get values for plot switches in focal plot routines
nodeflag=get(handles.plot_nodes_yes,'Value');
axflag=get(handles.plot_stress_yes,'Value');
dataflag=get(handles.plot_data_yes,'Value');

%  decide on dislocation vs moment tensor source models

stype=get(handles.use_dislocation_params,'Value');
%
% -------------------------------------------------------------------------
%  load up dislocation parameters, if wanted
if stype == 1
    mom(1,1)=str2num(get(handles.strike_input,'String'));
    mom(1,2)=str2num(get(handles.dip_input,'String'));
    mom(1,3)=str2num(get(handles.rake_input,'String'));
else
    %  load up moment tensor parameters
    mom(1,1)=str2num(get(handles.m11_input,'String'));
    mom(1,2)=str2num(get(handles.m22_input,'String'));
    mom(1,3)=str2num(get(handles.m33_input,'String'));
    mom(1,4)=str2num(get(handles.m12_input,'String'));
    mom(1,5)=str2num(get(handles.m13_input,'String'));
    mom(1,6)=str2num(get(handles.m23_input,'String'));
end
%
% -------------------------------------------------------------------------
%   send message to message strip
set(handles.message_strip,'String','Calculating Nodal Surfaces Plots');
%
% -------------------------------------------------------------------------
axes(handles.pwavenodes);
cla
set(handles.pwavenodes,'NextPlot','add','XTick',[],'YTick',[],'ZTick',[]);
pfocsphere(nmech,mom,stype,vp,vs,dens,p_pol,nodeflag,axflag,dataflag, z_displ)
%
% -------------------------------------------------------------------------
axes(handles.svwavenodes);
cla
set(handles.svwavenodes,'NextPlot','add','XTick',[],'YTick',[],'ZTick',[]);
svfocsphere(nmech,mom,stype,vp,vs,dens,sv_pol,nodeflag,axflag,dataflag, r_displ)
%
% -------------------------------------------------------------------------
axes(handles.shwavenodes);
cla
set(handles.shwavenodes,'NextPlot','add','XTick',[],'YTick',[],'ZTick',[]);
shfocsphere(nmech,mom,stype,vp,vs,dens,sh_pol,nodeflag,axflag,dataflag, t_displ)
%
% -------------------------------------------------------------------------
set(handles.pwavenodes,'NextPlot','replace');
set(handles.svwavenodes,'NextPlot','replace');
set(handles.shwavenodes,'NextPlot','replace');
%
% -------------------------------------------------------------------------
fig13 = figure (13);
clf 
   pfocsphere(nmech,mom,stype,vp,vs,dens,p_pol,nodeflag,axflag,dataflag, z_displ) 
   axis off;
   print(fig13, 'fig_P_pol.png', '-dpng')
   
% -------------------------------------------------------------------------
fig14 = figure (14);
clf 
   svfocsphere(nmech,mom,stype,vp,vs,dens,sv_pol,nodeflag,axflag,dataflag, r_displ)
   axis off;
   print(fig14, 'fig_SV_pol.png', '-dpng')
% -------------------------------------------------------------------------
fig15 = figure (15);
clf 
   shfocsphere(nmech,mom,stype,vp,vs,dens,sh_pol,nodeflag,axflag,dataflag, t_displ)
   axis off;
   print(fig15, 'fig_SH_pol.png', '-dpng')
% -------------------------------------------------------------------------
% Plotting the Swave polarization vector and P-nodal planes of the 
% best FMs on focalsphere
figure (16) 
clf 
    FM(1,1) = mom(1,1); FM(2,1) = mom(1,1);
    FM(1,2) = mom(1,2); FM(2,2) = mom(1,2);
    FM(1,3) = mom(1,3); FM(2,3) = mom(1,3);
    FM(1,4) = mom(1,3); FM(2,4) = mom(1,3)-1;
    
    plot_focal_figures(FM, 'sh_panel',p_pol,sv_pol,...
        sh_pol, vp, vs, dens,z_displ, r_displ, t_displ)
%
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------

% --- Executes on button press in pwave_plot_only.
function pwave_plot_only_Callback(hObject, eventdata, handles)
% hObject    handle to pwave_plot_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%   This function makes a separate P wave nodal surface plot
%
% ------------------------------------------------------------------------- 
%  set global variable for arrivals file
global vp vs dens nodeflag axflag dataflag stype mom
global p_pol nmech
%
% -------------------------------------------------------------------------
%   send message to message strip
set(handles.message_strip,'String','Calculating Separate P Nodal Plot');
%
% -------------------------------------------------------------------------
hp1=figure;
set(gca,'NextPlot','add','XTick',[],'YTick',[],'ZTick',[]);
pfocsphere(nmech,mom,stype,vp,vs,dens,p_pol,nodeflag,axflag,dataflag)
%
% -------------------------------------------------------------------------

% --- Executes on button press in svwave_plot_only.
function svwave_plot_only_Callback(hObject, eventdata, handles)
% hObject    handle to svwave_plot_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%   This function makes a separate SV wave nodal surface plot
% 
% -------------------------------------------------------------------------
%  set global variable for arrivals file
global vp vs dens nodeflag axflag dataflag stype mom
global sv_pol nmech
%
% -------------------------------------------------------------------------
%   send message to message strip
set(handles.message_strip,'String','Calculating Separate SV Nodal Plot');
%
% -------------------------------------------------------------------------
hsv1=figure;
set(gca,'NextPlot','add','XTick',[],'YTick',[],'ZTick',[]);
svfocsphere(nmech,mom,stype,vp,vs,dens,sv_pol,nodeflag,axflag,dataflag)
%
% -------------------------------------------------------------------------


% --- Executes on button press in shwave_plot_only.
function shwave_plot_only_Callback(hObject, eventdata, handles)
% hObject    handle to shwave_plot_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%   This function makes a separate SH wave nodal surface plot
%  
%  set global variable for arrivals file
%
global vp vs dens nodeflag axflag dataflag stype mom
global sh_pol nmech
%
% -------------------------------------------------------------------------
%   send message to message strip
set(handles.message_strip,'String','Calculating Separate SH Nodal Plot');
%
% -------------------------------------------------------------------------
hp1=figure;
set(gca,'NextPlot','add','XTick',[],'YTick',[],'ZTick',[]);
shfocsphere(nmech,mom,stype,vp,vs,dens,sh_pol,nodeflag,axflag,dataflag)
%
% -------------------------------------------------------------------------



























% Do not edit below this line

% --- Outputs from this function are returned to the command line.
function varargout = focal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function strike_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strike_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function strike_input_Callback(hObject, eventdata, handles)
% hObject    handle to strike_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strike_input as text
%        str2double(get(hObject,'String')) returns contents of strike_input as a double


% --- Executes during object creation, after setting all properties.
function dip_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dip_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dip_input_Callback(hObject, eventdata, handles)
% hObject    handle to dip_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dip_input as text
%        str2double(get(hObject,'String')) returns contents of dip_input as a double


% --- Executes during object creation, after setting all properties.
function rake_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rake_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function rake_input_Callback(hObject, eventdata, handles)
% hObject    handle to rake_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rake_input as text
%        str2double(get(hObject,'String')) returns contents of rake_input as a double


% --- Executes during object creation, after setting all properties.
function m11_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m11_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function m11_input_Callback(hObject, eventdata, handles)
% hObject    handle to m11_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m11_input as text
%        str2double(get(hObject,'String')) returns contents of m11_input as a double


% --- Executes during object creation, after setting all properties.
function m23_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m23_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function m23_input_Callback(hObject, eventdata, handles)
% hObject    handle to m23_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m23_input as text
%        str2double(get(hObject,'String')) returns contents of m23_input as a double


% --- Executes during object creation, after setting all properties.
function m13_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m13_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function m13_input_Callback(hObject, eventdata, handles)
% hObject    handle to m13_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m13_input as text
%        str2double(get(hObject,'String')) returns contents of m13_input as a double


% --- Executes during object creation, after setting all properties.
function m12_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m12_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function m12_input_Callback(hObject, eventdata, handles)
% hObject    handle to m12_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m12_input as text
%        str2double(get(hObject,'String')) returns contents of m12_input as a double


% --- Executes during object creation, after setting all properties.
function m33_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m33_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function m33_input_Callback(hObject, eventdata, handles)
% hObject    handle to m33_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m33_input as text
%        str2double(get(hObject,'String')) returns contents of m33_input as a double


% --- Executes during object creation, after setting all properties.
function m22_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m22_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function m22_input_Callback(hObject, eventdata, handles)
% hObject    handle to m22_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m22_input as text
%        str2double(get(hObject,'String')) returns contents of m22_input as a double


% --- Executes during object creation, after setting all properties.
function vp_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vp_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function vp_input_Callback(hObject, eventdata, handles)
% hObject    handle to vp_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vp_input as text
%        str2double(get(hObject,'String')) returns contents of vp_input as a double


% --- Executes during object creation, after setting all properties.
function vs_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vs_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function vs_input_Callback(hObject, eventdata, handles)
% hObject    handle to vs_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vs_input as text
%        str2double(get(hObject,'String')) returns contents of vs_input as a double


% --- Executes during object creation, after setting all properties.
function dens_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dens_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dens_input_Callback(hObject, eventdata, handles)
% hObject    handle to dens_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dens_input as text
%        str2double(get(hObject,'String')) returns contents of dens_input as a double


% --- Executes on button press in plot_data_yes.
function plot_data_yes_Callback(hObject, eventdata, handles)
% hObject    handle to plot_data_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_data_yes

set(handles.plot_data_yes,'Value',1);
set(handles.plot_data_no,'Value',0);


% --- Executes on button press in plot_data_no.
function plot_data_no_Callback(hObject, eventdata, handles)
% hObject    handle to plot_data_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_data_no

set(handles.plot_data_yes,'Value',0);
set(handles.plot_data_no,'Value',1);


% --- Executes on button press in plot_nodes_yes.
function plot_nodes_yes_Callback(hObject, eventdata, handles)
% hObject    handle to plot_nodes_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_nodes_yes

set(handles.plot_nodes_yes,'Value',1);
set(handles.plot_nodes_no,'Value',0);


% --- Executes on button press in plot_nodes_no.
function plot_nodes_no_Callback(hObject, eventdata, handles)
% hObject    handle to plot_nodes_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_nodes_no

set(handles.plot_nodes_yes,'Value',0);
set(handles.plot_nodes_no,'Value',1);


% --- Executes on button press in plot_stress_yes.
function plot_stress_yes_Callback(hObject, eventdata, handles)
% hObject    handle to plot_stress_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_stress_yes

set(handles.plot_stress_yes,'Value',1);
set(handles.plot_stress_no,'Value',0);

% --- Executes on button press in plot_stress_no.
function plot_stress_no_Callback(hObject, eventdata, handles)
% hObject    handle to plot_stress_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_stress_no

set(handles.plot_stress_yes,'Value',0);
set(handles.plot_stress_no,'Value',1);



% --- Executes during object creation, after setting all properties.
function message_strip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to message_strip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function message_strip_Callback(hObject, eventdata, handles)
% hObject    handle to message_strip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of message_strip as text
%        str2double(get(hObject,'String')) returns contents of message_strip as a double


% --- Executes on button press in use_dislocation_params.
function use_dislocation_params_Callback(hObject, eventdata, handles)
% hObject    handle to use_dislocation_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_dislocation_params

set(handles.use_dislocation_params,'Value',1);
set(handles.use_moment_params,'Value',0);

% --- Executes on button press in use_moment_params.
function use_moment_params_Callback(hObject, eventdata, handles)
% hObject    handle to use_moment_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_moment_params

set(handles.use_dislocation_params,'Value',0);
set(handles.use_moment_params,'Value',1);



% --- Executes on button press in choose_arrivals.
function choose_arrivals_Callback(hObject, eventdata, handles)
% hObject    handle to choose_arrivals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function strike_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strike_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function strike_inc_Callback(hObject, eventdata, handles)
% hObject    handle to strike_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strike_inc as text
%        str2double(get(hObject,'String')) returns contents of strike_inc as a double


% --- Executes during object creation, after setting all properties.
function dip_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dip_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dip_inc_Callback(hObject, eventdata, handles)
% hObject    handle to dip_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dip_inc as text
%        str2double(get(hObject,'String')) returns contents of dip_inc as a double


% --- Executes during object creation, after setting all properties.
function rake_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rake_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function rake_inc_Callback(hObject, eventdata, handles)
% hObject    handle to rake_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rake_inc as text
%        str2double(get(hObject,'String')) returns contents of rake_inc as a double


% --- Executes during object creation, after setting all properties.
function number_mech_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_mech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function number_mech_Callback(hObject, eventdata, handles)
% hObject    handle to number_mech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_mech as text
%        str2double(get(hObject,'String')) returns contents of number_mech as a double






% % p_amp = [6.00804e3 6.02727e2 7.3e2 3.6023e2 1.091e2 -2.1e2];
% % 1.74726 -4.73478 -2.06493e1];
% % sv_amp = [-3.93701e3 2.61682e3 -1.42244e4 -2.5429e3 1.56539e2 1.1022e3];
% % sh_amp = [1.27738e4 4.55512e3 -8.80765e3 -1.816e3 -5.0783e2 -1.86125e3];
% 
% p_amp = [6.00804e3 6.02727e2 7.3e2 3.6023e2 1.091e2 -2.1e2 1.74726 -4.73478 -2.06493e1];
% weight_P_amp = [5 5 5 5 5 5 5 5 5];
% % 1.74726 -4.73478 -2.06493e1];
% sv_amp = [-7934.81626 5635.33484 -12820.86034 -2276.4942 62 -808.09277 0 0 0];
% sh_amp = [1.27811e4 4.56889e3 -9.1082e3 3.00613e3 -488 -3.7086e3 0 0 0 ];
% weight_S_amp = [5 5 5 5 5 5 0 0 0];




% path_to_SAC = '/usr/local/sac/bin';
% path_to_travelt = '/Users/oluwaseunfadugba/Documents/Waveform_Modeling_Now/time';
% % -------------------------------------------------------------------------
% % Inputs (Real data)
% sta_no = [1        2        3        4         5         6         7        8        9];
% dist   = [9.424156 25.1711  35.6124  44.78913  65.8215   68.58533  331.5408 332.98   367.2564];
% az     = [23.3895  126.803  73.69663 38.52594  40.56081  57.86335  241.3568 90.34291 60.07621];
% 
% % Syntax: data = [station_no p_amp weight];
% % p_data = [1     6.00804e3   5;
% %           2     6.02727e2   5;
% %           3     7.3e2       5;
% %           4     3.6023e2    5;
% %           5     1.091e2     5;
% %           6     -2.1e2      5];
% 
% p_data = [1     6.00804e3   5;
%           2     6.02727e2   5;
%           3     7.3e2       5;
%           4     3.6023e2    5;
%           5     1.091e2     5;
%           6     -2.1e2      5;
%           7     1.74726     5;
%           8     -4.73478    5;
%           9     -2.06493e1  5];
% 
% sv_data = [1    -7934.81626     5;
%            2    5635.33484      5;
%            3    -12820.86034    5;
%            4    -2276.4942      5;
%            5    62              5;
%            6    -808.09277      5];
% 
% sh_data = [1    1.27811e4       5;
%            2    4.56889e3       5;
%            3    -9.1082e3       5;
%            4    3.00613e3       5;
%            5    -488            5;
%            6    -3.7086e3       5];
% 
% 
% use_sv_pol = 0; % 0 -- do not use sv polarity
%                 % 1 -- use sv polarity
% use_sv_amp = 3; % use_sv_amp = 0 -- do not use sv at all
%                 % use_sv_amp = 1 -- use only |sv|/sh
%                 % use_sv_amp = 2 -- use only |sv|/p
%                 % use_sv_amp = 3 -- use both |sv|/sh and |sv|/p
% 
% w_pol = 0.5; w_ratio = 0.5;
% 
% depth = 24.5; depth_error = 4;
% halfspace = 0; synthetic = 1; 
% 
% % syntax: vel_model =[vp vs rho depth layer_no]
% vel_model =[6.08   3.51   2.732    0.0    1;
%             6.25   3.60   2.756    6.0   1; 
%             6.55   3.70   2.789    12.0   1;
%             7.2    4.20   2.90    40.0   1;
%             8.0    4.60   3.2    40.0001   2;
%             8.0    4.60   3.2    200.0   2]; 
%     
% %  vel_model = [6.0    3.5   2.7    0.0    1; ... % Halfspace model
% %               6.0    3.5   2.7    40.0   1];
% %
