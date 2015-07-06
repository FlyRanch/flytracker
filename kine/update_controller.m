% Function UPDATE_CONTROLLER
% 
% CALLING FUNCTION: load_images, get_calibration
% ACTIONS: update controller display to reflect images and/or calibration
%          loaded
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function update_controller

controller = findobj('Tag','controller');
data = guidata(controller);

% Clear all objects with tag "obj_control"
to_delete = findobj('Tag','obj_control');
delete(to_delete)

% If both images and calibration are loaded also...
if isfield(data,'images') & isfield(data,'cal')
    if isfield(data,'config') == 0
        object_dialog                               % Get user to choose object configuration (sets kinematics uicontrols)
    end
    
    show_objects('get_kinematics')                  % Make get kinematics uicontrols visible
    show_objects('obj_controls')                    % Make frame for object controls visible
    toggle_pointer_fcn('Mark Points')               % Set pointer function to mark points
else
    set(data.handles.get_kinematics,'Visible','off')
    set(data.handles.obj_controls,'Visible','off')
end

% If the images are loaded
if isfield(data,'images')
    set_text('image_text','[''Current image sequence: '',data.images.file]',[])     % Display path/file names in text strings
    show_objects('slider')
else
    set_text('image_text','< No Image Sequence Loaded >')
    set_text('data_text','< Data Not Saved >')
    set(data.handles.slider,'Visible','off') 
end

% If the calibration is loaded
if isfield(data,'cal')
    if isfield(data.cal,'filename')                                      % to grandfather in old saved data
        set_text('cal_text',['Calibration loaded:',data.cal.filename])   % Display Calibration text
    else
        set_text('cal_text','Calibration successful')
    end
else
    set_text('cal_text','< No Calibration Loaded >') 
end

