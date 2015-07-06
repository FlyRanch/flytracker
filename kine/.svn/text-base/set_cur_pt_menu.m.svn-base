% Function SET_CUR_PT_MENU
% 
% CALLING FUNCTION: Callback for cur_object_menu in controller;
%                   set_cur_obj_menu
% ACTIONS: Sets cur_point_menu based on fields in data.kine.(object)
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 26, 2004 by gwyneth

function set_cur_pt_menu

controller = findobj('Tag','controller');
data = guidata(controller);

% What object is selected?
object = get_string(data.handles.cur_object_menu);

% Get points for current object
points = data.kine.(object).config.points;
if isempty(points)
    points = {'NONE'};
end

% Set cur_point_menu
set(data.handles.cur_point_menu,'String',points)
set(data.handles.cur_point_menu,'Value',1)%%%%DEBUG

% Set other menus dependent on current object
set_obj_param_menus

% Clear all the object controls
to_delete = findobj('Tag','obj_control');
delete(to_delete)

% Set the environment to be correct for the current object type
type = data.kine.(object).config.type;
mfile_path = data.setup.mfile_path;
cd([mfile_path,filesep,'object_types',filesep,type])
dig_setup(object)
cd(mfile_path)


% Change coordinate plot to be correct point
cur_pt_menu_callback
%     plot_it  %%%DEBUG - need to find right order to do things
%     draw_it
