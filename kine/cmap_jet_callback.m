% Function CMAP_JET_CALLBACK
% 
% CALLING FUNCTION: Callback for cmap jet button
% ACTIONS: resets colormap for dig_fig
% PARENT PROGRAM: Kine_v2_1
% LAST MODIFIED: April 29, 2005 by astraw

function cmap_jet_callback

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

axes(dig_data.handles.(['cam',num2str(1)]))
colormap('jet')