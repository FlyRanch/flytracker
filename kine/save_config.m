% Function SAVE_CONFIG
% 
% CALLING FUNCTION: Callback for save_button in obj_fig
% ACTIONS: Saves current configuration (all objects) to a .mat file in
%          Kine_v2_0\object_configs, reloads the config menu, sets value to
%          recently saved config
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 21, 2004 by gwyneth

function save_config

obj_data = guidata(gcf);

choice = get(obj_data.handles.config_menu,'Value');
num_configs = length(obj_data.configs);

if choice > num_configs
    current_config = 'custom';
else
    current_config = obj_data.configs{choice};
end

objects = fieldnames(obj_data.kine);
vars = [];

if strcmp(current_config,'custom') == 0
    warndlg('That configuration is already saved','Saved Configuration')
    return
end

if length(objects) == 0
    warndlg('There are no objects to save','No Objects')
    return
end

config_name = inputdlg('New configuration name:','Save Custom Configuration');
if isempty(config_name) | strcmp(config_name,'')
    warndlg('No configuration name entered.  No configuration was saved.','No Configuration Saved')
    return
end

save_path = [obj_data.config_path,config_name{1},'.mat'];

for i = 1:length(objects)
    % Make the variables to be saved (make variables with object names outside
    % kine structure) - want only the config field saved
    command = [objects{i},'.config = obj_data.kine.(objects{i}).config;'];
    eval(command)
    
    % Make string of names of variables to be saved (starts with comma)
    vars = [vars,',''',objects{i},''''];
end

already_there = strmatch(config_name,obj_data.configs);
while already_there ~= 0
    overwrite = questdlg('This configuration name is already being used, do you want to overwrite it?','Save Custom Configuration');
    switch overwrite
        case 'Cancel'
            return
            
        case 'No'
            config_name = inputdlg('New configuration name:','Save Custom Configuration');
            if isempty(config_name) | strcmp(config_name,'')
                warndlg('No configuration name entered.  No configuration was saved.','No Configuration Saved')
                return
            end
            
            save_path = [obj_data.config_path,config_name{1},'.mat'];
            already_there = strmatch(config_name,obj_data.configs);
            
        case 'Yes'
            break
    end
end

command = ['save(''',save_path,'''',vars,')'];
eval(command)

config_menu_CreateFcn
obj_data = guidata(gcf);

val = strmatch(config_name,obj_data.configs);
set(obj_data.handles.config_menu,'Value',val)
% load_object_config                        % This just confirms that the
% configuration is saved appropriately -- have commented out so that won't
% overwrite current coordinate data when save config 5/12/04
