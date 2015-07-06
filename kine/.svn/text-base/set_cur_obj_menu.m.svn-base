% Function SET_CUR_OBJ_MENU
% 
% CALLING FUNCTION: obj_ok_callback, load_data
% ACTIONS: Sets cur_object_menu and calls function to set cur_point_menu 
%          in controller based on fields in data.kine
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 26, 2004 by gwyneth

function set_cur_obj_menu

controller = findobj('Tag','controller');
data = guidata(controller);

objects = fieldnames(data.kine);
set(data.handles.cur_object_menu,'String',objects)
set(data.handles.cur_object_menu,'Value',1)%%%%%%%%5DEBUG

set_cur_pt_menu
% set_obj_param_menus
