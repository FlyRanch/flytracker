% Function COLOR_LIST_CALLBACK
% 
% CALLING FUNCTION: Callback for color_list listbox
% ACTIONS: Sets data.kine.(cur_obj).color equal to selected string
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 3, 2004 by gwyneth

function color_list_callback

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);   % Get current object
old_color = data.kine.(cur_obj).config.color;                % Get current type
new_color = get_string(data.handles.color_list);      % Get selected type

if strcmp(old_color,new_color)                        % Do nothing if same color clicked on
    return
end

data.kine.(cur_obj).config.color = new_color;
guidata(controller,data)

draw_it
calc_it
plot_it
