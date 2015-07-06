% Function GET_LISTBOX_STRING
% 
% CALLING FUNCTION: edit_object, delete_object  
% ACTIONS: In Object Dialog: object listbox items are specially formatted strings combining
%          several variables about an object.  These strings are made in
%          set_listbox_string.  This function "knows" the format used to make them
%          and can extract any information required from them.  Should only
%          need to get object name - then can find other info in
%          obj_data.kine.(obj_name).
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function obj_name = get_listbox_string

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

% set_listbox_string format:
%   listbox_string{i} = [objects{i},'  (',type,', ',num_pts,' points, ',color,', visible: ',vis,')'];

% Get name of object selected in obj_fig listbox
listed_objects = get(obj_data.handles.obj_list,'String');
choice = get(obj_data.handles.obj_list,'Value');
obj_string = listed_objects{choice};

obj_name = '';
for i = 1:length(obj_string)
    s = obj_string(i);
    if strcmp(s,' ')==0
        obj_name = [obj_name,s];
    else
        break
    end
end
        
