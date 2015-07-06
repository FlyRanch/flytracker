% Function sel_string = GET_STRING(menu_handle)
% 
% CALLING FUNCTION: type_list_callback
% ACTIONS: Get string from selected choice in popupmenu or listbox
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 3, 2004 by gwyneth

function sel_string = get_string(menu_handle)

options = get(menu_handle,'String');
choice = get(menu_handle,'Value');
sel_string = options{choice};