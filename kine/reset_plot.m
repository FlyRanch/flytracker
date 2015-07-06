% Function RESET_PLOT
% 
% CALLING FUNCTION: callback for reset_plot button in dig_fig
% ACTIONS: autoscale axes, replot red triangle
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 30, 2004 by gwyneth

function reset_plot

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

axes(dig_data.handles.coord_plot)
axis auto
plot_it

% Can't find the bug, so query radio buttons, set pointer function to agree
% with them
% switch get(dig_data.handles.goto_radio,'Value')
%     
%     case get(dig_data.handles.goto_radio,'Max')
%         toggle_pointer_fcn('Mark Points')
%         
%     case get(dig_data.handles.goto_radio,'Min')
%         toggle_pointer_fcn('Zoom')
%         
% end
