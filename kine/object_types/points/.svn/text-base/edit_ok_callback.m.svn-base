% Function EDIT_OK_CALLBACK
% 
% CALLING FUNCTION: Kine_v2_0
% ACTIONS: Edit dialog okay button callback for type 'POINTS'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 7, 2004 by gwyneth

function edit_ok_callback

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

obj_name = get(edit_data.handles.name_edit,'String');
pt_edit_handles = findobj('Tag','pt_edit');
% NTOE: Want pt_edit_handles in order, so will have to sort them by position of
% the edit boxes with LARGEST y values first (position(2))

point_names = {};

for i = 1:length(pt_edit_handles)
    pos = get(pt_edit_handles(i),'Position');
    point_order(i) = pos(2);
    
    point_names{i} = get(pt_edit_handles(i),'String');
    
    if isempty(point_names{i})
        warndlg('There was an unnamed point!')
    end
    
end

[sorted, ind] = sort(point_order);   % sorted from low y to high y
point_names(ind) = point_names;
point_names = fliplr(point_names);  % so need to reverse order

obj_data.kine.(obj_name).config.points = point_names;
guidata(obj_fig,obj_data)
