% Function KEYBOARD_COMMANDS
% 
% CALLING FUNCTION: KeyPressFcn for dig_fig windows 
% ACTIONS: Performs the correct action depending on which key pressed
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 25, 2004 by gwyneth

function keyboard_commands

key = get(gcf,'CurrentCharacter');

% don't allow any keyboard commands if in the middle of digitizing a point:
dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

if length(dig_data.temp) > 1
    return
end

switch key
    
    case 'r' % Next Point
        advance
        
    case 'e' % Previous Point
        advance
        
    case '`' % Esc mark point (clear lines)
        % Delete all objects in dig_data.temp
        points_to_delete = fieldnames(dig_data.temp);
        for i = 1:length(points_to_delete)
            delete(dig_data.temp.(points_to_delete{i}))
            dig_data.temp = rmfield(dig_data.temp,points_to_delete{i});
        end
        
        guidata(dig_fig,dig_data)
        toggle_pointer_fcn('Mark Points')
        
    case 'p' % Clear current point
        clear_pt('p')
        
    case 'o' % Clear current object
        clear_pt('o')
        
    case 'f' % Next frame
        h = findobj('Tag','frame_box');
        frame = eval(get(h,'String'));
        new_frame = frame + 1;
        
        set(h,'String',num2str(new_frame))
        
        update_images

    case 'd' % Previous frame
        h = findobj('Tag','frame_box');
        frame = eval(get(h,'String'));
        new_frame = frame - 1;
        
        set(h,'String',num2str(new_frame))
        
        update_images
        
    case 'g' % Next frame
        h = findobj('Tag','frame_box');
        frame = eval(get(h,'String'));
        new_frame = frame + 5;
        
        set(h,'String',num2str(new_frame))
        
        update_images

    case 's' % Previous frame
        h = findobj('Tag','frame_box');
        frame = eval(get(h,'String'));
        new_frame = frame - 5;
        
        set(h,'String',num2str(new_frame))
        
        update_images
        
    case 'x' % Set pointer function to Mark Points
        toggle_pointer_fcn('Mark Points')
        
    case 'z' % Set pointer function to Zoom
        toggle_pointer_fcn('Zoom')
        
    case '=' % Save
        save_data
        
    case 'm' % Open object dialog
        open_object_dialog
        
    case 'q' % Reset images (unzoom)
        reset_button_callback
        
    case 'a' % Load images
        load_images
        
    case ';' % Calibrate
        load_calibration
        
    case 'l' % Load Data
        load_data
        
    case 'v' % Toggle current object visibility
        h = findobj('Tag','vis_check');
        cur_val = get(h,'Value');
        if cur_val == get(h,'Min')
            new_val = get(h,'Max');
        else
            new_val = get(h,'Min');
        end
        set(h,'Value',new_val)
        
        vis_check_callback
        
    case 'b' % Toggle autosave on/off
        h = findobj('Tag','autosave_check');
        cur_val = get(h,'Value');
        if cur_val == get(h,'Min')
            new_val = get(h,'Max');
        else
            new_val = get(h,'Min');
        end
        set(h,'Value',new_val)

        save_auto
        
        
    case 'y'
        inc_amount = 0.001;
        inc_dec(inc_amount,'x')
        
    case 'h'
        inc_amount = -0.001;
        inc_dec(inc_amount,'x')
        
    case 'u'
        inc_amount = 0.001;
        inc_dec(inc_amount,'y')
        
    case 'j'
        inc_amount = -0.001;
        inc_dec(inc_amount,'y')
        
    case 'i'
        inc_amount = 0.001;
        inc_dec(inc_amount,'z')
        
    case 'k'
        inc_amount = -0.001;
        inc_dec(inc_amount,'z')
        
end



%===========================================================
%COMMENT ON THE FOLLOWING FUNCTIONS!!!!
%===========================================================
function inc_dec(inc_amount, curraxis)

controller = findobj('Tag','controller');
data = guidata(controller);

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

switch(curraxis)
    case 'x'
        axis = 1;
    case 'y'
        axis = 2;
    case 'z'
        axis = 3;
    
    % defaults to x-axis
    otherwise
        curraxis = 'x';
        axis = 1;
end

coords(axis) = coords(axis) + inc_amount;
    
store(coords)
draw_it
calc_it
plot_it
