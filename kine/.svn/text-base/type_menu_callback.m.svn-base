% Function TYPE_MENU_CALLBACK
% 
% CALLING FUNCTION: Type menu in edit object dialog
% ACTIONS: Looks at type menu to see which type selected, then runs the
%          type_menu_callback script inside that type's folder
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 4, 2004 by gwyneth

function type_menu_callback

controller = findobj('Tag','controller');
data = guidata(controller);

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

mfile_path = data.setup.mfile_path;
type = get_string(edit_data.handles.type_menu);

set(edit_data.handles.num_edit,'Style','edit')

run([mfile_path,filesep,'object_types',filesep,type,filesep,'type_menu_callback'])