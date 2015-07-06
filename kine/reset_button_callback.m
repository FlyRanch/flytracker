% Function RESET_BUTTON_CALLBACK
% 
% CALLING FUNCTION: Callback for reset button
% ACTIONS: resets zoom for dig_fig axes
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function reset_button_callback

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

for i = 1:data.setup.cam_num
    axes(dig_data.handles.(['cam',num2str(i)]))
    zoom out
end