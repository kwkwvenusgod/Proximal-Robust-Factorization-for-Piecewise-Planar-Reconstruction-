function varargout = view_outputs(varargin)
% VIEW_OUTPUTS MATLAB code for view_outputs.fig
%      VIEW_OUTPUTS, by itself, creates a new VIEW_OUTPUTS or raises the existing
%      singleton*.
%
%      H = VIEW_OUTPUTS returns the handle to a new VIEW_OUTPUTS or the handle to
%      the existing singleton*.
%
%      VIEW_OUTPUTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_OUTPUTS.M with the given input arguments.
%
%      VIEW_OUTPUTS('Property','Value',...) creates a new VIEW_OUTPUTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_outputs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_outputs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose 'GUI allows only one
%      instance to run (singleton)'.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_outputs

% Last Modified by GUIDE v2.5 20-Sep-2012 17:37:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @view_outputs_OpeningFcn, ...
    'gui_OutputFcn',  @view_outputs_OutputFcn, ...
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


% --- Executes just before view_outputs is made visible.
function view_outputs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure

% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_outputs (see VARARGIN)

% Choose default command line output for view_outputs
handles.output = hObject;


% to modify
handles.outpath = get(handles.editRoot,'String');
handles.outpath ='~/';
if handles.outpath(end) ~='/'
    handles.outpath= [handles.outpath '/'];
end
%{
video_names={'camera_motion1','car1_stabilized','car2_stabilized',...
    'dog','phone','phone','phone','sample1','sample3','sample4',...
    'sample5','birdfall2','cheetah','girl','monkeydog','parachute',...
    'penguin','monkeydog'};
video_id = [1:length(video_names)];
handles.video_id = containers.Map(video_names, video_id);
%}

handles.video_names = get(handles.lstFiles,'String');
handles.vid = get(handles.lstFiles,'Value');
handles.video_name = handles.video_names{handles.vid};

handles.btn_tags = {'btn_kpos' 'btn_kapp' 'btn_area' 'btn_alpha' 'btn_beta','btn_deltap_scale'};
for i =1:numel(handles.btn_tags)
    if eval(['get(handles.' handles.btn_tags{i} ',''Value'')'])==1
        handles.btn_id = i;
        break;
    end
end



handles.kpos = [100,1000,10000,100000];
handles.kapp = [100,1000,10000,100000];
handles.area = [400];
handles.alpha = [-15,-10,-5];
handles.beta = [-3,-2,-1];
handles.deltap_scale = [0.25,0.5,1,2];
handles.para_size = [length(handles.kpos),length(handles.kapp),...
    length(handles.area),length(handles.alpha),length(handles.beta),length(handles.deltap_scale)];






% precompute all stats variable
if ~exist('output_stats.mat')
    num_para = prod(handles.para_size);
    num_video = numel(handles.video_names);
    stats=zeros(num_video,num_para,8);
    
    for i = 1: num_para
        % loop through all parameters
        for j = 1: num_video
            % loop through all videos
            try
                tmp = load([handles.outpath num2str(i,'%06d') '/seg_' handles.video_names{j}],'stats');
                stats(i,j,:) = tmp.stats;
            end
        end
    end
    stats(:,:,1) = stats(:,:,1)/1000;
    stats(:,:,4) = stats(:,:,4)/100;
    stats(:,:,7) = stats(:,:,7)/50;
    stats(:,:,8) = stats(:,:,8)/500;
    save output_stats.mat stats
else
    load('output_stats.mat','stats');
end
handles.stats = stats;
delete stats
%scale the data so that they can be plotted in one figure

handles.stats_tag={'sv num (x1000)','ue 2D','accu 2D','br 2D (x100)','ue 3D','accu 3D','br 3D (x50)','H (x500)'};
handles.stats_color={'r.-','g.-','b.-','c.-','y.-','cx-','rx-','gx-'};

% parameter range correspondence
handles.stats_ran = cell(1,numel(handles.para_size));
handles.p_id = ones(1,numel(handles.para_size));
handles.output_id = '000001';


handles=Set_SPlabel(handles);
handles=Set_Img(handles);
Set_SlideBar(handles);
handles.img_id=round(get(handles.slider_img,'Value'));
Set_PlotSP(handles)
Set_PlotStats(handles);

% configurations
try
    set(handles.slider_kpos, 'Value', 1);
    set(handles.slider_kpos, 'Max', handles.para_size(1));
    set(handles.slider_kpos, 'Min', 1);
    set(handles.slider_kpos, 'SliderStep',[1,1]/(handles.para_size(1)-1));
end
try
    set(handles.slider_kapp, 'Value', 1);
    set(handles.slider_kapp, 'Max', handles.para_size(2));
    set(handles.slider_kapp, 'Min', 1);
    set(handles.slider_kapp, 'SliderStep',[1,1]/(handles.para_size(2)-1));
end
try
    set(handles.slider_area, 'Value', 1);
    set(handles.slider_area, 'Max', handles.para_size(3));
    set(handles.slider_area, 'Min', 1);
    set(handles.slider_area, 'SliderStep',[1,1]/(handles.para_size(3)-1));
end
try
    set(handles.slider_alpha, 'Value', 1);
    set(handles.slider_alpha, 'Max', handles.para_size(4));
    set(handles.slider_alpha, 'Min', 1);
    set(handles.slider_alpha, 'SliderStep',[1,1]/(handles.para_size(4)-1));
end
try
    set(handles.slider_beta, 'Value', 1);
    set(handles.slider_beta, 'Max', handles.para_size(5));
    set(handles.slider_beta, 'Min', 1);
    set(handles.slider_beta, 'SliderStep',[1,1]/(handles.para_size(5)-1));
end

try
    set(handles.slider_deltap_scale, 'Value', 1);
    set(handles.slider_deltap_scale, 'Max', handles.para_size(6));
    set(handles.slider_deltap_scale, 'Min', 1);
    set(handles.slider_deltap_scale, 'SliderStep',[1,1]/(handles.para_size(6)-1));
end

% Update handles structure
guidata(hObject, handles);









% UIWAIT makes view_outputs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view_outputs_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function handles = Set_SPlabel(handles)
handles=Parse_para2id(handles);
%img configuration
%sp configuration
try
    tmp = load([handles.outpath handles.output_id '/seg_' handles.video_name],'sp_labels');
    handles.sp_labels= tmp.sp_labels;
    handles.sp_len=size(handles.sp_labels,3);
catch
    handles.sp_labels=[];
    handles.sp_len=2;
    disp('NO SP_label available')
end

function handles=Set_Img(handles)
try
    tmp = load([handles.outpath , handles.output_id '/seg_' handles.video_name],'folder');
    %handles.imgfolder = tmp.folder;
    handles.imgfolder = ['/home/Stephen/Desktop/Bird/Segment/data' tmp.folder(find(tmp.folder=='t',1,'first')+6:end)];
    handles.imgs = U_getims(handles.imgfolder);
catch
    handles.imgs=[];
end


function handles=Set_SlideBar(handles)
set(handles.slider_img, 'Value', 1);
set(handles.slider_img, 'Max', handles.sp_len);
set(handles.slider_img, 'Min', 1);
set(handles.slider_img, 'SliderStep', [1,1]/(handles.sp_len-1));





function handles=Parse_para2id(handles)
% plot stat ids
for i =1:6
    new_ind=num2cell(handles.p_id);
    new_ind{7-i}=1:handles.para_size(7-i);
    handles.stats_ran{i} = ...
        (new_ind{6}-1)*prod(handles.para_size(1:end-1))+...
        (new_ind{5}-1)*prod(handles.para_size(1:end-2))+...
        (new_ind{4}-1)*prod(handles.para_size(1:end-3))+...
        (new_ind{3}-1)*prod(handles.para_size(1:end-4))+...
        (new_ind{2}-1)*prod(handles.para_size(1:end-5))+...
        new_ind{1};
end

% parameter id
tmp_id = handles.stats_ran{handles.btn_id}(handles.p_id(1));
handles.output_id = num2str(tmp_id, '%06d');
handles.output_id
%{
% test parameter id extraction
[cov_var_p, cov_var_a, area_var, alpha, beta_min_alpha] = ndgrid(...
    [100,1000,10000,100000],...
    [100,1000,10000,100000],...
    [400],...
    [-1000:100:-400],...
    [100:100:900]);
[cov_var_p(para_id); cov_var_a(para_id);...
area_var(para_id); alpha(para_id); beta_min_alpha(para_id)]
%}


function Set_PlotSP(handles)
axes(handles.axes1);
hold off
try
    imagesc(handles.sp_labels(:,:,handles.img_id))
    %{
    transparent =0.5;
    im = handles.imgs{handles.img_id};
    N = size(im,1)*size(im,2);
    mask = find(ismember(handles.sp_labels(:,:,handles.img_id), handles.sp_len));
    im(mask) = im(mask)*transparent;
    im(mask+N) = im(mask+N)*transparent + (1-transparent);
    im(mask+2*N) = im(mask+2*N)*transparent;
    borders = is_border_valsIMPORT(double(handles.sp_labels(:,:,handles.img_id)));
    im = setPixelColors(im, find(borders), [1 1 0]);
    imagesc(im)
    %}
catch
    text(0,0,'No Image available')
    %imagesc(handles.sp_labels(:,:,handles.img_id))
end
axis off;



function Set_PlotStats(handles)
axes(handles.axes2);
%need to scale data properly

% replot it every time since not expensive...
% (local update of the current parameter may need some work to make it look nice)
cla
hold on
para_id = handles.stats_ran{handles.btn_id};
for i = 1:length(handles.stats_tag)
    plot(1:numel(para_id),handles.stats(handles.vid,para_id,i),handles.stats_color{i})
end
legend(handles.stats_tag)

% highlight current parameter set
try
    x= find(para_id==str2num(handles.output_id));
    plot(x, squeeze(handles.stats(handles.vid,x,:)),'ko')
catch
    disp('Can''t find current parameter settings')
end


axis on






% --- Executes on selection change in lstFiles.
function lstFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lstFiles (see GCBO)

% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstFiles
handles.vid = get(handles.lstFiles,'Value');
if ~strcmp(handles.video_name,handles.video_names{handles.vid})
    handles.video_name = handles.video_names{handles.vid};
    handles=Set_SPlabel(handles);
    handles=Set_Img(handles);
    Set_SlideBar(handles);
    Set_PlotSP(handles)
    Set_PlotStats(handles)
end
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_kpos_Callback(hObject, eventdata, handles)
handles.p_id(1) = round(get(hObject,'Value'));
set(handles.text_kpos,'String',num2str(handles.kpos(handles.p_id(1))))
handles = Set_SPlabel(handles);
Set_PlotSP(handles)
Set_PlotStats(handles)
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_kapp_Callback(hObject, eventdata, handles)
handles.p_id(2) = round(get(hObject,'Value'));
set(handles.text_kapp,'String',num2str(handles.kapp(handles.p_id(2))))
handles = Set_SPlabel(handles);
Set_PlotSP(handles);
Set_PlotStats(handles)
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_area_Callback(hObject, eventdata, handles)
handles.p_id(3) = round(get(hObject,'Value'));
set(handles.text_area,'String',num2str(handles.area(handles.p_id(3))))
handles=Set_SPlabel(handles);
Set_PlotSP(handles);
Set_PlotStats(handles)
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_alpha_Callback(hObject, eventdata, handles)
handles.p_id(4) = round(get(hObject,'Value'));
set(handles.text_alpha,'String',num2str(handles.alpha(handles.p_id(4))))
handles=Set_SPlabel(handles);
Set_PlotSP(handles);
Set_PlotStats(handles)
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_beta_Callback(hObject, eventdata, handles)
handles.p_id(5) = round(get(hObject,'Value'));
set(handles.text_beta,'String',num2str(handles.beta(handles.p_id(5))))
handles=Set_SPlabel(handles);
Set_PlotSP(handles);
Set_PlotStats(handles)
guidata(hObject, handles);

function slider_deltap_scale_Callback(hObject, eventdata, handles)
handles.p_id(6) = round(get(hObject,'Value'));
set(handles.text_deltap_scale,'String',num2str(handles.deltap_scale(handles.p_id(6))))
handles=Set_SPlabel(handles);
Set_PlotSP(handles);
Set_PlotStats(handles)
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_img_Callback(hObject, eventdata, handles)
handles.img_id = round(get(hObject,'Value'));
Set_PlotSP(handles)
Set_PlotStats(handles)
guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
if ~isempty(eventdata.OldValue),          % Check for an old selected object
    oldTag = get(eventdata.OldValue,'Tag');   % Get Tag of old selected object
    index = find(strcmp(oldTag,handles.btn_tags));     % Find index of match in buttonTags
    switch index
        case 1
            set(handles.btn_kpos,'Value',0);
        case 2
            set(handles.btn_kapp,'Value',0);
        case 3
            set(handles.btn_area,'Value',0);
        case 4
            set(handles.btn_alpha,'Value',0);
        case 5
            set(handles.btn_beta,'Value',0);
        case 6
            set(handles.btn_deltap_scale,'Value',0);
    end
end

newTag = get(eventdata.NewValue,'Tag');   % Get Tag of new selected object
handles.btn_id = find(strcmp(newTag,handles.btn_tags));     % Find index of match in buttonTags
switch handles.btn_id
    case 1
        set(handles.btn_kpos,'Value',1);
    case 2
        set(handles.btn_kapp,'Value',1);
    case 3
        set(handles.btn_area,'Value',1);
    case 4
        set(handles.btn_alpha,'Value',1);
    case 5
        set(handles.btn_beta,'Value',1);
    case 6
        set(handles.btn_deltap_scale,'Value',1);
end

Set_PlotStats(handles)
guidata(hObject,handles);  % Save handles structure


function editRoot_Callback(hObject, eventdata, handles)
handles.outpath = get(handles.editRoot,'String');
if handles.outpath(end) ~='/'
    handles.outpath= [handles.outpath '/'];
end
guidata(hObject,handles);  % Save handles structure


























function lstFiles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_kpos_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_kapp_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_area_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_alpha_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_beta_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_deltap_scale_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_img_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function editRoot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text1_CreateFcn(hObject, eventdata, handles)
function text2_CreateFcn(hObject, eventdata, handles)
function text3_CreateFcn(hObject, eventdata, handles)
function text4_CreateFcn(hObject, eventdata, handles)
function text5_CreateFcn(hObject, eventdata, handles)
function text6_CreateFcn(hObject, eventdata, handles)
function text_kpos_CreateFcn(hObject, eventdata, handles)
function text_kapp_CreateFcn(hObject, eventdata, handles)
function text_area_CreateFcn(hObject, eventdata, handles)
function text_alpha_CreateFcn(hObject, eventdata, handles)
function text_beta_CreateFcn(hObject, eventdata, handles)
function text_deltap_scale_CreateFcn(hObject, eventdata, handles)
