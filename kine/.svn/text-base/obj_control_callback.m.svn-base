% Function OBJ_CONTROL_CALLBACK(fcn)
% 
% CALLING FUNCTION: callback for object-specific uicontrols in the
%                   controller
% ACTIONS: calls correct callback function from correct type folder;
%          assumes that the name of the callback function is stored as a string in a
%          cell array in the objects 'UserData'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: October 2, 2004 by gwyneth


function obj_control_callback(fcn)

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);
cur_type = data.kine.(cur_obj).config.type;

mfile_path = data.setup.mfile_path;
run([mfile_path,filesep,'object_types',filesep,cur_type,filesep,fcn])