% Function TYPE_MENU_CALLBACK
% 
% CALLING FUNCTION: Kine_v2_0
% ACTIONS: Type menu callback for type 'POINTS'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 7, 2004 by gwyneth

function type_menu_callback

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

% Clear all objects specific to type, resize, set num_edit to 0   
to_clear = findobj('Parent',edit_fig);

for i = 1:length(edit_data.handles.to_keep)
    ind = find(to_clear ~= edit_data.handles.to_keep(i));
    to_clear = to_clear(ind);
end

delete(to_clear)

pos = get(edit_data.handles.num_edit,'Position');
y_sub = pos(2) - 6;

fig_pos = get(edit_data.handles.edit_fig,'Position');
fig_pos(2) = fig_pos(2) + y_sub; % move down by y_add
fig_pos(4) = fig_pos(4) - y_sub; % make longer by y_add
set(edit_data.handles.edit_fig,'Position',fig_pos)

for j = 1:length(edit_data.handles.to_move)
    pos = get(edit_data.handles.to_move(j),'Position');
    pos(2) = pos(2) - y_sub; % move down by y_add
    set(edit_data.handles.to_move(j),'Position',pos)
end

set(edit_data.handles.num_edit,'String','0')
set(edit_data.handles.num_edit,'Callback','num_edit_callback')
