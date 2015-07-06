% Function SHOW_OBJECTS(obj_handles)
% 
% CALLING FUNCTION: load_images, get_calibration
% ACTIONS: Makes visible the objects with handles listed in data.handles.obj_handles.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 10, 2004 by gwyneth
       
function show_objects(obj_handles)

controller = findobj('Tag','controller');
data = guidata(controller);

set(data.handles.(obj_handles),'Visible','on')
