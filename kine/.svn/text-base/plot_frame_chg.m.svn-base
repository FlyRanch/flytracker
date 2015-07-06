% Function PLOT_FRAME_CHG
% 
% CALLING FUNCTION: open_dig_fig, plot_it
% ACTIONS: When the user clicks on the coordinate plot, changes the current
%          frame to be that corresponding to the x-value of the clicked point
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 30, 2004 by gwyneth

function plot_frame_chg

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

cur_pt = get(dig_data.handles.coord_plot,'CurrentPoint');

new_frame = cur_pt(1); % See Axes Properties, CurrentPoint for the output matrix of CurrentPoint
new_frame = round(new_frame);

set(data.handles.frame_box,'String',new_frame)
update_images