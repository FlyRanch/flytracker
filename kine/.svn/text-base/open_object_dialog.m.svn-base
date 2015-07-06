% Function OPEN_OBJECT_DIALOG
% 
% CALLING FUNCTION: Callback for edit_obj_button in controller
% ACTIONS: calls object_dialog then sets custom parameters if data.config
%          is set to 'custom'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 26, 2004 by gwyneth

function open_object_dialog

controller = findobj('Tag','controller');
data = guidata(controller);

object_dialog

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

if strcmp(data.config,'custom')
    add_custom
else
    config = strmatch(data.config,obj_data.configs);
    set(obj_data.handles.config_menu,'Value',config)
    obj_data.cur_config = config;
end

obj_data.kine = data.kine;
guidata(obj_fig,obj_data)

set_listbox_string
