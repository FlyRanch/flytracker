% Function LOAD_DATA
% 
% CALLING FUNCTION: load_data button callback
% ACTIONS: Loads data for a certain trial and the images associated with it
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 1, 2004 by gwyneth

function load_data

controller = findobj('Tag','controller');
cur_data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
if isempty(dig_fig) == 0    % if there is already a dig_fig window
    dig_fig_closereq        % then close it
    
    dig_fig = findobj('Tag','dig_fig');
    if isempty(dig_fig) == 0    % if there is still a dig_fig window then the user cancelled the action
        return
    end
end

uiload

% check loaded data to make sure has all necessary fields
if exist('data') == 0 | isfield(data,'setup') == 0 |...
        isfield(data,'save') == 0 | isfield(data,'advance') == 0 |...
        isfield(data,'colors') == 0 | isfield(data,'images') == 0 |...
        isfield(data.cal.coeff,'DLT_1') == 0 | isfield(data,'kine') == 0 |...
        isfield(data,'config') == 0
    warndlg('The selected file did not contain a recognizeable data structure','No data found')
    return
end
        
dig_fig = findobj('Tag','dig_fig');
if isempty(dig_fig) == 0
    dig_fig_closereq
end

%--------------------------------------------------------------------------
% FOR COMPATIBILITY with data saved with old versions of kine that had
% different setup fields, check to see that the default setup file has all 
% the correct fields (from Kine_v2_1); adjust paths that may have changed
setup_wrong = 0;
if isfield(data.setup,'name') == 0
    data.setup.name = 'ByHand';
end
if isfield(data.setup,'cam_num') == 0
    ans = inputdlg('How many cameras do you have data from?')
    data.setup.cam_num = str2num(ans{1});
end
if isfield(data.setup,'cal_type') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.cal_type = 'DLT';
end
if isfield(data.setup,'cal_mfile') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.cal_mfile = 'get_calibration_points_g';
end
if isfield(data.setup,'image_filter') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.image_filter = '*.bmp';
end
if isfield(data.setup,'mfile_path') == 0
    data.setup.mfile_path = uigetdir('','Locate the Kine folder');
elseif exist(data.setup.mfile_path) == 0
    disp('The location of your Kine folder has changed since the data were last saved')
    data.setup.mfile_path = uigetdir('','Locate the Kine folder');
end
if isfield(data.setup,'data_path') == 0
    data.setup.data_path = uigetdir('','Locate the data folder');
elseif exist(data.setup.data_path) == 0
    disp('The location of data folder no longer exists')
    data.setup.data_path = uigetdir('','Locate the data folder');
    
    disp('Click on any image from the desired video sequence')
    [temp, data.images.path] = uigetfile('*.bmp');
end
if isfield(data.setup,'keyboard_command_file') == 0
    ans = inputdlg('What is the name of your keyboard command file?');
    data.setup.keyboard_command_file = ans{1};
end

if setup_wrong ~= 0
    msg = ['Hard-wired setup values were used for ',num2str(setup_wrong),' parameter(s).'];
    dlgname = 'HARD-WIRED SETUP VALUES USED';
    warndlg(msg,dlgname)
end

%--------------------------------------------------------------------------
% FOR COMPATIBILITY with old way of loading, find image file root here if
% doesn't already exist (this is code from get_image_file_info.m), 20040820

if isfield(data.images,'file_root') == 0
    namesStruct = dir( [data.images.path , data.setup.image_filter] );    % struct array :)
    namesStruct_length = length( namesStruct );                 % number of files in each folder
        filenameCells = {namesStruct(:).name};
    filenameCells = sort( filenameCells );
    root = 1;
    n = 0;
    while root == 1
        n = n + 1;
        root = strncmp(filenameCells{1}, filenameCells{end}, n);
    end
    data.images.file_root = data.images.file(1:n-1);
end
%--------------------------------------------------------------------------

data.handles = cur_data.handles;    % get new handles structure
guidata(controller,data)            % set data to newly loaded data

load_images         % load corresponding images
set_cur_obj_menu    % set menus
save_auto           % Set autosave appropriately
%draw_it             % Draw stored kine points on current frame

show_objects('get_kinematics')                  % Make get kinematics uicontrols visible
set_text('cal_text','Calibration successful')   % Display Calibration text
set_text('data_text',['Last save to ''',data.save.filename,''' at ',data.save.timestamp])

