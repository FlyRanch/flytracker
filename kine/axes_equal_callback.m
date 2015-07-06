% Function axes_equal_CALLBACK
% 
% CALLING FUNCTION: Callback for reset button
% ACTIONS: sets image aspect ratio to 1-to-1 for dig_fig axes
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 4, 2004 by doug

function axes_equal_callback

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

for i = 1:data.setup.cam_num
    axes(dig_data.handles.(['cam',num2str(i)]))
    axis image
end