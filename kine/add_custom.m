% Function ADD_CUSTOM
% 
% CALLING FUNCTION: add_object, edit_ok_callback, delete_object, make_all_vis,
%                   order_ok_callback
% ACTIONS: Adds 'custom' string to the config_menu if obj_fig, not saved
%          out to obj_data, though, so when change config_menu selection (run
%          config_menu callback) the 'custom' option goes away.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 21, 2004 by gwyneth

function add_custom

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

obj_data.cur_config = length(obj_data.configs) + 1;
guidata(obj_fig,obj_data)

% Update config_menu to include "custom" and make this the current setting
I = strmatch('custom',obj_data.configs,'exact');
if isempty(I)
    obj_data.configs{end+1} = 'custom';
end

set(obj_data.handles.config_menu,'String',obj_data.configs,'Value',length(obj_data.configs))
