% Function CLEAR_PT(char)
% 
% CALLING FUNCTION: Called by keyboard commands, 'o' and 'p'
% ACTIONS: Deletes either the current point or the current object
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 26, 2004 by gwyneth

function clear_pt(char)

controller = findobj('Tag','controller');
data = guidata(controller);

frame = str2num(get(data.handles.frame_box,'String'));
cur_pt = get_string(data.handles.cur_point_menu);
cur_obj = get_string(data.handles.cur_object_menu);

switch char
    
    case 'p'
        proceed = questdlg(['Are you sure you want to clear the point ''',cur_pt,'?'],'Clear Object','Yes','Cancel','Yes');
        switch proceed
            case 'Yes'
                % Only clear current frame!!
                num_pt = get(data.handles.cur_point_menu,'Value');           % get the number of the point
                data.kine.(cur_obj).data.coords(num_pt,:,frame) = zeros(1,3);    % set point values to zero
                guidata(controller,data)
            case 'Cancel'
                toggle_pointer_fcn('Mark Points')
                return
        end
        
    case 'o'
        proceed = questdlg(['Are you sure you want to clear the object ''',cur_obj,'?'],'Clear Object','Yes','Cancel','Yes');
        switch proceed
            case 'Yes'
                % Only clear current frame!!
                num_pts = length(data.kine.(cur_obj).config.points);               % get the number of points for this object
                data.kine.(cur_obj).data.coords(:,:,frame) = zeros(num_pts,3);   % make a zeros array for the first frame (3 for 3 dimensions)
                guidata(controller,data)
            case 'Cancel'
                toggle_pointer_fcn('Mark Points')
                return
        end
        
end

draw_it
calc_it
plot_it