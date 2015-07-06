% Function OBJ_OK_CALLBACK
% 
% CALLING FUNCTION: Callback for ok_button in obj_fig
% ACTIONS: Saves obj_data.kine structure into data.kine structure according
%          to how much has been changed (may delete already digitized points)
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 26, 2004 by gwyneth

function obj_ok_callback

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

controller = findobj('Tag','controller');
data = guidata(controller);

% Determine which configuration is chosen
choice = get(obj_data.handles.config_menu,'Value');
if choice > length(obj_data.configs)
    new_config = 'custom';
else
    new_config = obj_data.configs{choice};
end

% Get cell array of all object names
new_objects = fieldnames(obj_data.kine);

% Check that the array is not empty otherwise don't allow
if isempty(new_objects)
    warndlg('Can''t have a configuration with no objects.','No Objects')
    return
end

% Make data.kine and data.config
data.kine = obj_data.kine;
data.config = new_config;

guidata(controller,data)
delete(obj_data.handles.obj_fig)

% Figure out which objects are new, run dig_setup for these (should be all
% initially)
new_ones = intersect(new_objects,fieldnames(data.kine)); % data.kine might not exist yet...

for i = 1:length(new_ones)
    % Set the environment to be correct for the current object type
    type = obj_data.kine.(new_ones{i}).config.type;
    mfile_path = data.setup.mfile_path;
    cd([mfile_path,filesep,'object_types',filesep,type])
    dig_setup(new_ones{i})
    cd(mfile_path)
end

set_cur_obj_menu
