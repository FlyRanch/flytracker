% Function EDIT_OBJECT
% 
% CALLING FUNCTION: Callback for Edit button in obj_fig for 'ARC'
%                   objects
% ACTIONS: Opens edit object dialog then sets values of uicontrols
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 10, 2004 by gwyneth

function edit_object

addpath C:\MATLAB6p5\work\Kine_v2_0  %%%DEBUG!!!

% Set value of uicontrols based on currently selected object
obj_fig = findobj('Tag','obj_fig');
edit_fig = findobj('Tag','edit_fig');

obj_data = guidata(obj_fig);
edit_data = guidata(edit_fig);

obj_name = get(edit_data.handles.name_edit,'String');

% Set num_edit callback correctly for arc
type_menu_callback

% Set number of points
num_pts = obj_data.kine.(obj_name).num_pts;
set(edit_data.handles.num_edit,'String',num2str(num_pts))

