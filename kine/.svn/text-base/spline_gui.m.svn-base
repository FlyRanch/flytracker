function varargout = spline_gui(varargin)
% SPLINE_GUI M-file for spline_gui.fig
%      SPLINE_GUI, by itself, creates a new SPLINE_GUI or raises the existing
%      singleton*.
%
%      H = SPLINE_GUI returns the handle to a new SPLINE_GUI or the handle to
%      the existing singleton*.
%
%      SPLINE_GUI('Property','Value',...) creates a new SPLINE_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to spline_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SPLINE_GUI('CALLBACK') and SPLINE_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SPLINE_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spline_gui

% Last Modified by GUIDE v2.5 28-Jun-2004 05:13:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spline_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @spline_gui_OutputFcn, ...
                   'gui_LayoutFcn',  @spline_gui_LayoutFcn, ...
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


% --- Executes just before spline_gui is made visible.
function spline_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for spline_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spline_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spline_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function txtStartFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function txtStartFrame_Callback(hObject, eventdata, handles)
% hObject    handle to txtStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtStartFrame as text
%        str2double(get(hObject,'String')) returns contents of txtStartFrame as a double


% --- Executes during object creation, after setting all properties.
function txtEndFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function txtEndFrame_Callback(hObject, eventdata, handles)
% hObject    handle to txtEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEndFrame as text
%        str2double(get(hObject,'String')) returns contents of txtEndFrame as a double


% --- Executes on button press in btnOK.
function btnOK_Callback(hObject, eventdata, handles)
% hObject    handle to btnOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

start_frame = str2num(get(handles.txtStartFrame, 'String'));
end_frame = str2num(get(handles.txtEndFrame, 'String'));

if (start_frame < 1)
    msgbox('Start frame must be greater than 0','Error','error')
elseif (end_frame <= start_frame)
    msgbox('End frame must be greater than start frame','Error','error')
elseif (end_frame > data.images.frames)
    msgbox(['End frame must be less than total number of frames (' num2str(data.images.frames) ')'],'Error','error')
else
    spline_data(start_frame, end_frame);
    delete(gcf)
end


% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)



% --- Creates and returns a handle to the GUI figure. 
function h1 = spline_gui_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

h1 = figure(...
'Units','characters',...
'Color',[0.937254901960784 0.937254901960784 0.937254901960784],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','spline_gui',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[103.8 40.6923076923077 75.2 20.7692307692308],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',zeros(1,0));

setappdata(h1, 'GUIDEOptions', struct(...
'active_h', 1.020033e+002, ...
'taginfo', struct(...
'figure', 2, ...
'edit', 3, ...
'text', 4, ...
'pushbutton', 4), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'lastSavedFile', 'C:\Documents and Settings\Qing\Desktop\spline_gui.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','spline_gui(''txtStartFrame_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[9.8 7.53846153846154 20.2 1.61538461538462],...
'String','',...
'Style','edit',...
'CreateFcn','spline_gui(''txtStartFrame_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','txtStartFrame');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','spline_gui(''txtEndFrame_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[39.8 7.53846153846154 20.2 1.61538461538462],...
'String','',...
'Style','edit',...
'CreateFcn','spline_gui(''txtEndFrame_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','txtEndFrame');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[9.8 9.07692307692308 20.2 1.53846153846154],...
'String','Start frame',...
'Style','text',...
'Tag','lblStartFrame');


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[39.8 9.07692307692308 20.2 1.53846153846154],...
'String','End frame',...
'Style','text',...
'Tag','lblEndFrame');


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontName','Tahoma',...
'FontSize',12,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[9.8 12.9230769230769 50.2 2.76923076923077],...
'String','Spline Tool',...
'Style','text',...
'Tag','text3');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','spline_gui(''btnOK_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[27.6 3.07692307692308 15 2.23076923076923],...
'String','OK',...
'Tag','btnOK');


h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','spline_gui(''btnCancel_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[45 3.07692307692308 15 2.23076923076923],...
'String','Cancel',...
'Tag','btnCancel');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      SPLINE_GUI, by itself, creates a new SPLINE_GUI or raises the existing
%      singleton*.
%
%      H = SPLINE_GUI returns the handle to a new SPLINE_GUI or the handle to
%      the existing singleton*.
%
%      SPLINE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPLINE_GUI.M with the given input arguments.
%
%      SPLINE_GUI('Property','Value',...) creates a new SPLINE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.4 $ $Date: 2002/05/31 21:44:31 $

gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [getfield(gui_State, gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % SPLINE_GUI
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % SPLINE_GUI('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % SPLINE_GUI(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA
        guidata(gui_hFigure, guihandles(gui_hFigure));
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) & ischar(varargin{ind+1}) & ...
                strncmpi(varargin{ind},'visible',len1) & len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index
            break;
        end
        try, set(gui_hFigure, varargin{index}, varargin{index+1}), catch, break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)
if nargin('openfig') == 3 
    gui_hFigure = openfig(name, singleton, 'auto');
else
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end

