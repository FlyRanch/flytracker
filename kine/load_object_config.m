% Function LOAD_OBJECT_CONFIG
% 
% CALLING FUNCTION: callback for config_menu in object dialog; save_config;
%                   createFcn for obj_list listbox in object dialog
% ACTIONS: loads appropriate .mat file from object_config folder, sets
%          object dialog listbox appropriately, when "okay" clicked this will set
%          the data.kine structures appropriately
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 18, 2004 by gwyneth

function load_object_config

controller = findobj('Tag','controller');
data = guidata(controller);

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

config_menu = findobj('Tag','config_menu'); % Need to find the object because 
                                            % guihandles not yet created when this 
                                            % called as the config_menu creation function
choice = get(config_menu,'Value');
current_config = obj_data.configs{choice};

% In the case that load_object_config is being called as the callback for
% the config_menu, if a config has previously been loaded, check to see
% that the user wants to change config and clear any existing data
if isfield(data,'config') & gcbo == findobj('Tag','config_menu')
    proceed = questdlg('Changing the configuration will clear the data.  Do you want to proceed?','Change Configuration','Yes','No','Yes');
    if strcmp(proceed,'No')
        set(config_menu,'Value',obj_data.cur_config)
        return
    end
end

kine = load([obj_data.config_path,current_config,'.mat']);

objects = fieldnames(kine);

if isempty(objects)
    obj_data.kine = struct([]);
else
    obj_data.kine = kine;
end

obj_data.cur_config = choice;
guidata(gcf,obj_data)

% Reset config menu options to get rid of custom option if it's there
set(config_menu,'String',obj_data.configs,'Value',choice)
set_listbox_string
