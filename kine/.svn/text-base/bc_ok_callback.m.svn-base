% Function BC_OK_CALLBACK
% 
% CALLING FUNCTION: Callback for OK button edit in bc_fig
% ACTIONS: Reads the values of all the uicontrols and outputs these to the
%          .bc field of the data structure
% PARENT PROGRAM: Kine_v2_1
% LAST MODIFIED: November 21, 2006 by gwyneth

function bc_ok_callback

bc_data = guidata(gcf);

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

bc_fig = findobj('Tag','bc_fig');
bc_data = guidata(bc_fig);

data.bc.bc_obj = get_string(bc_data.handles.bc_obj_menu);

data.bc.axis1.object1 = get_string(bc_data.handles.axis1_obj1_menu);
data.bc.axis1.point1 = get_string(bc_data.handles.axis1_point1_menu);
data.bc.axis1.object2 = get_string(bc_data.handles.axis1_obj2_menu);
data.bc.axis1.point2 = get_string(bc_data.handles.axis1_point2_menu);
data.bc.axis1.length = [];

data.bc.axis2.object1 = get_string(bc_data.handles.axis2_obj1_menu);
data.bc.axis2.point1 = get_string(bc_data.handles.axis2_point1_menu);
data.bc.axis2.object2 = get_string(bc_data.handles.axis2_obj2_menu);
data.bc.axis2.point2 = get_string(bc_data.handles.axis2_point2_menu);
data.bc.axis2.length = [];

data.bc.origin.object = get_string(bc_data.handles.axis3_obj3_menu);
data.bc.origin.point = get_string(bc_data.handles.axis3_point3_menu);

guidata(controller, data)

% Make BC views visible
set(dig_data.handles.bc_view,'Visible','on')

% Delete bc figure
delete(gcf)
