% Function VIS_CHECK_CALLBACK
% 
% CALLING FUNCTION: Callback for vis_check checkbox in controller
% ACTIONS: Sets 'visible' field of current object to 'on' or 'off'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 3, 2004 by gwyneth

function vis_check_callback

controller = findobj('Tag','controller');
data = guidata(controller);

checked = get(data.handles.vis_check,'Value');

switch checked
    
    case get(data.handles.vis_check,'Max') % box checked
        vis = 'on';
        
    case get(data.handles.vis_check,'Min') % box not checked
        vis = 'off';
        
end

cur_obj = get_string(data.handles.cur_object_menu); % get current object
data.kine.(cur_obj).config.visible = vis;

guidata(controller,data)

draw_it