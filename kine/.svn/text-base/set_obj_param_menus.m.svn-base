% Function SET_OBJ_PARAM_MENUS
% 
% CALLING FUNCTION: obj_ok_callback, load_data, plot_it
% ACTIONS: Sets cur_object_menu and calls function to set cur_point_menu 
%          in controller based on fields in data.kine
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 26, 2004 by gwyneth

function set_obj_param_menus

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);

% Type listbox
val = strmatch(data.kine.(cur_obj).config.type, get(data.handles.type_list,'String'));
    if isempty(val)
        val = [];
    end
set(data.handles.type_list,'Value',val)
    
% Color listbox
val = strmatch(data.kine.(cur_obj).config.color, get(data.handles.color_list,'String'));
    if isempty(val)
        val = [];
    end
set(data.handles.color_list,'Value',val)

% Visible checkbox
switch data.kine.(cur_obj).config.visible
    case 'on'
        set(data.handles.vis_check,'Value',get(data.handles.vis_check,'Max'))
    case 'off'
        set(data.handles.vis_check,'Value',get(data.handles.vis_check,'Min'))
end

