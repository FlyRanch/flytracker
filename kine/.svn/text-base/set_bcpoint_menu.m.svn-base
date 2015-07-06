% Function SET_BEPOINT_MENU
% 
% CALLING FUNCTION: Callback for object menus bc_fig dialog
% ACTIONS: Sets point_menu based on fields in data.kine.(object)
% PARENT PROGRAM: Kine_v2_
% LAST MODIFIED: November 21, 2006 by gwyneth

function set_bcpoint_menu

controller = findobj('Tag','controller');
data = guidata(controller);

bc_fig = findobj('Tag','bc_fig');
bc_data = guidata(bc_fig);

obj_name = get(gco,'Tag');
points = data.kine.(get_string(bc_data.handles.(obj_name))).config.points;

ax_num = obj_name(5);
obj_num = obj_name(10);

point_name = ['axis',ax_num,'_point',obj_num,'_menu'];

set(bc_data.handles.(point_name),'String',points)