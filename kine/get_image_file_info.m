 % Function success = GET_IMAGE_FILE_INFO
% 
% CALLING FUNCTION: load_images
% ACTIONS: gets image path and file name inputs from user
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 10, 2004 by gwyneth

function success = get_image_file_info
controller = findobj('Tag','controller');
data = guidata(controller);    % Retrieve guidata

% if gcbo == data.handles.load_data_button
if isfield(data,'images')
    %check to see if those images still exist in the right place
    if exist([data.images.path,data.images.file]) ~= 0
        success = 1;
        return
    else
        disp('The images corresponding to the requested data could not be found in the location they were last saved.  Please locate the correct image files.');
    end
end

% Save path and file names in guidata structure
cd(data.setup.data_path)

[ data.images.file , data.images.path ] = uigetfile(data.setup.image_filter,'Select any frame from desired video sequence'); 
    if data.images.file == 0
        success = 0;
        cd(data.setup.mfile_path)
        return;
    else
        success = 1;
    end
    
cd(data.setup.mfile_path)

% Auto-detect the number of cameras with data available
dir_str = findstr(filesep,data.images.path);     
exp_path = data.images.path(1:dir_str(end-1));  % get path one folder up

cam_list = dir( [ exp_path, 'cam*' ] );         % assume cam folders start with "cam"
cams_found = length(cam_list);

if cams_found == 0
    % Save path and file names in guidata structure
    data.images.flexload(1).file = data.images.file;
    data.images.flexload(1).path = data.images.path;

    for i = 2:data.setup.cam_num
        cd(data.setup.data_path)

        [ data.images.flexload(i).file , data.images.flexload(i).path ] = uigetfile(data.setup.image_filter,['Select any frame from Camera ',num2str(i)]);
        if data.images.flexload(i).file == 0
            success = 0;
            cd(data.setup.mfile_path)
            return;
        else
            success = 1;
        end

        cd(data.setup.mfile_path)
    end

    for i = 1:data.setup.cam_num
        filenameCells = makeFileList( data.images.flexload(i).path, data.setup.image_filter );
        data.images.flexload(i).file_root = findRoot( filenameCells );
        data.images.frames = findFrames( filenameCells, data.images.flexload(i).file_root );
    end
    
else

    if cams_found > data.setup.cam_num
        question = ['Kine is currently setup to display ',...
            num2str(data.setup.cam_num),...
            ' cameras, but you have data for ',...
            num2str(cams_found),...
            ' cameras.  Would you like to change your camera setup to display all the data?'];

        change_cam_num = questdlg(question, 'Camera Setup', 'YES','NO','YES');

    elseif cams_found < data.setup.cam_num
        change_cam_num = 'YES';

    elseif cams_found == data.setup.cam_num
        change_cam_num = 'NO';

    end

    switch change_cam_num
        case 'YES'
            data.setup.cam_num = cams_found;
            guidata(controller,data)
        case 'NO' % no change needed
    end
    
    filenameCells = makeFileList( data.images.path, data.setup.image_filter );
    data.images.file_root = findRoot( filenameCells );
    data.images.frames = findFrames( filenameCells, data.images.file_root );
    
end

guidata(controller,data)   % Store guidata     

%--------------------------------------------------------------------------
function filenameCells = makeFileList( path, image_filter )

% Get list of all .bmp files in directory
namesStruct = dir( [path , image_filter] );    % struct array :)
namesStruct_length = length( namesStruct );                 % number of files in each folder

% Sort the list
filenameCells = {namesStruct(:).name};
filenameCells = sort( filenameCells );

%--------------------------------------------------------------------------
function file_root = findRoot( filenameCells )

% Find the root of the file name (without frame values appended) - the only
% assumption made is that the name is of the form [root, numbers].bmp
% Names should be in order, so compare first and last to get value range
% All letters 1:n-1 are the root, n:end are the frame numbers

root = 1;
n = 0;
while root == 1

    n = n + 1;
    root = strncmp(filenameCells{1}, filenameCells{end}, n);

end

file_root = filenameCells{1}(1:n-1);

%--------------------------------------------------------------------------
function frames = findFrames( filenameCells, file_root );

% HARDWIRED: length of appended file type is 4
ft_length = 4;
n = length(filenameCells{1}) - length(file_root) - ft_length;

% Assume last four string values are the suffix (usually .bmp)
first_frame = str2num(filenameCells{1}(end-(ft_length+n-1):end-ft_length));
last_frame = str2num(filenameCells{end}(end-(ft_length+n-1):end-ft_length));

frames = last_frame - first_frame + 1;
