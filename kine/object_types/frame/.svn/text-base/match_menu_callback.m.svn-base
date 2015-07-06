% Function MATCH_MENU_CALLBACK
% 
% CALLING FUNCTION: callback for match object menu in edit_fig
% ACTIONS: sets appropriate match point menu with correct choices
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: October 2, 2004 by gwyneth

function match_menu_callback

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

% Find out which menu was used
menu_tag = get(gcbo,'Tag');

% Get corresponding point menu handle
h_pt_menu = findobj('Tag',['pt_menu_',menu_tag(end)]);

% Get the new point names for this menu
cur_obj = get_string(gcbo);
pts = obj_data.kine.(cur_obj).config.points;

% Set the point menu appropriately
set(h_pt_menu,'String',pts)