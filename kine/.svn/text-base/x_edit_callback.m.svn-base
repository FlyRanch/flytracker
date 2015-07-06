% Function:       x_edit_callback
% 
% Description:    This function is called when the edit box for the x axis
%                 is changed by the user and updates the point that is
%                 currently being edited.
%                 
% Arguments:      None.
%             
% Modifies:       data - structure containing all digitized data
% Returns:        None.
% Error handling: None.
%
% Revision history:   07/06/2004  Qing Liu - initial revision

function x_edit_callback

controller = findobj('Tag','controller');
data = guidata(controller);

edit_value = str2num(get(data.handles.x_edit,'String'));
set(data.handles.x_slider,'Value',edit_value);

% Get current object,  point, and frame from menus
choice = get(data.handles.cur_object_menu,'Value');
objects = get(data.handles.cur_object_menu,'String');
cur_object = objects{choice};                               % Note that while the cur_object is retrieved as a string
                                                            % because it will be used to retrieve a structure field, the
cur_point = get(data.handles.cur_point_menu,'Value');       % cur_point is retrieved as a number because it will be used
                                                            % to index the rows of a matrix - this number should be okay without
frame = str2num(get(data.handles.frame_box,'String'));      % verification because the menu is made by retrieving the points
                                                            % string, so the name of the point will be data.kine.(object).points{cur_point}
coords = data.kine.(cur_object).data.coords(cur_point,:,frame);

coords(1) = edit_value;
    
store(coords)

calc_it
plot_it
draw_it