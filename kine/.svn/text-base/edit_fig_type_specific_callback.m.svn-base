% Function EDIT_FIG_TYPE_SPECIFIC_CALLBACK
% 
% CALLING FUNCTION: callback for edit_fig uicontrols created by
%                   'type_menu_callback'
% ACTIONS: calls correct callback function from correct type folder;
%          assumes that the name of the callback function is stored as a string in a
%          cell array in the objects 'UserData'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: October 2, 2004 by gwyneth


function edit_fig_type_specific_callback

controller = findobj('Tag','controller');
data = guidata(controller);

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

cur_type = get_string(edit_data.handles.type_menu);
fcn = get(gcbo,'UserData');

mfile_path = data.setup.mfile_path;
run([mfile_path,filesep,'object_types',filesep,cur_type,filesep,fcn{1}])