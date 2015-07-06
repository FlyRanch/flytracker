% Function DISPLAY_IMAGE
% 
% CALLING FUNCTION: load_images
% ACTIONS: Reads and displays correct image frame for each camera
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 10, 2004 by gwyneth

function display_image

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

missed_images = 0;

% Get current path, file names and frame number from data structure
frame = get(data.handles.frame_box,'String');
file_type = data.setup.image_filter(2:end);%'.bmp';

% See if there is a data.images.flexload field, if so load these values in
% the camera loop
flexload = isfield(data.images,'flexload');

% Load new images into memory:
for i = 1:data.setup.cam_num

    if flexload
        path_name = data.images.flexload(i).path;
        file_name = data.images.flexload(i).file;
        root = data.images.flexload(i).file_root;
    else
        path_name = data.images.path;
        file_name = data.images.file;
        root = data.images.file_root;
        % For each camera, change path name to look in correct folder
        path_name(end-1) = num2str(i);
    end

    % For each frame, change file name to get correct file
    digits_to_add = length(file_name) - length(root) - length(file_type);
    num_zeros = digits_to_add - length(frame);
    
    zero_string = '';    % Make string of zeros to append to frame #
    for j = 1:num_zeros
        zero_string = [zero_string,'0'];
    end

    % OPTION: use code here if root is different by one digit for the 3
    % different cameras (hard-wired for Mike Tu's data for now)
    % root(4) = num2str(i);
    
    file_to_load = [ root, zero_string, frame, file_type ];
    
    % Set image name to read
    full_name = [ path_name, file_to_load];
    
    % Check to make sure the image exists
    if exist(full_name) == 0
        missed_images = missed_images +1;
        continue
    end

    % Read in image
    cam(i).theImage = imread( full_name,'bmp' );
    colormap gray
end

% Delete old images and display new ones
for i = 1:data.setup.cam_num
    
    % Find and delete old images
    old_images = findobj('Tag',['im',num2str(i)]);                                  
    %delete(old_images)
    if isempty(old_images) == 0
        %dig_data.handles = rmfield(dig_data.handles,['im',num2str(i)]);
        set(dig_data.handles.(['im',num2str(i)]),'CData',cam(i).theImage)
        set(dig_data.handles.(['cam',num2str(i)]),'Tag',['cam',num2str(i)],'FontSize',8,'NextPlot','add','DrawMode','fast')
    else
    % Display image    
    axes(dig_data.handles.(['cam',num2str(i)]));      % Set axes to appropriate camera (1,2,3)
    dig_data.handles.(['im',num2str(i)]) = imagesc(cam(i).theImage,'Tag',['im',num2str(i)],'EraseMode','none');   % Set image handle to 'im1', 'im2', 'im3'   
    set(gca,'Tag',['cam',num2str(i)],'FontSize',8,'NextPlot','add','DrawMode','fast','Box','off','TickLength',[0 0])%'XTick',[],'YTick',[]);
    end
    %delete(old_images) % moved deleting of old images down here (after displaying new ones)
                        % so don't get ugly flash on image change; didn't
                        % need to do this in matlab6.5 for some reason;
                        % image redraw is still slow; note this workds
                        % because handle of image to be deleted is taken
                        % from dig_data.handles and put into old_images

end

if missed_images > 0
    warndlg('Not all video frames can be found.  Check that images are stored with the same prefix under different folders for each camera');
end
    
guidata(dig_fig,dig_data); 