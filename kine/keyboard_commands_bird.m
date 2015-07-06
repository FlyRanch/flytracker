% Function KEYPRESS
% 
% CALLING FUNCTION: KeyPressFcn for dig_fig windows 
% ACTIONS: Performs the correct action depending on which key pressed
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 21, 2004 by Qing Liu


% WARNING: This assumes a particular configuration for points to be
% digitized (i.e. item 1 on object menu is the left wing, item 1 on point
% menu is hinge for wing objects, head for body objects).

function keyPress

controller = findobj('Tag','controller');
data = guidata(controller);

key = get(gcf,'CurrentCharacter');

switch key
    
    % set current point to left wingtip
    case 'a'
        set(data.handles.cur_object_menu, 'Value', 1)
        set_cur_pt_menu
        set(data.handles.cur_point_menu, 'Value', 2)
        toggle_pointer_fcn('Mark Points')      
      
    % set current point to left wing hinge
    case 's'
        set(data.handles.cur_object_menu, 'Value', 1)
        set_cur_pt_menu
        set(data.handles.cur_point_menu, 'Value', 1)
        toggle_pointer_fcn('Mark Points')      
        
    % set current point to right wing hinge
    case 'd'
        set(data.handles.cur_object_menu, 'Value', 3)
        set_cur_pt_menu
        set(data.handles.cur_point_menu, 'Value', 1)
        toggle_pointer_fcn('Mark Points')      
        
    % set current point to right wingtip
    case 'f'
        set(data.handles.cur_object_menu, 'Value', 3)
        set_cur_pt_menu
        set(data.handles.cur_point_menu, 'Value', 2)
        toggle_pointer_fcn('Mark Points')      
        
    % set current point to head
    case 'e'
        set(data.handles.cur_object_menu, 'Value', 2)
        set_cur_pt_menu
        set(data.handles.cur_point_menu, 'Value', 1)
        toggle_pointer_fcn('Mark Points')      
        
    % set current point to tail
    case 'c'
        set(data.handles.cur_object_menu, 'Value', 2)
        set_cur_pt_menu
        set(data.handles.cur_point_menu, 'Value', 2)
        toggle_pointer_fcn('Mark Points')      
        

    case 'm'
        dig_fig = findobj('Tag','dig_fig');
        dig_data = guidata(dig_fig);

        image_axes = [];
		for i = 1:data.setup.cam_num
            image_axes = [image_axes dig_data.handles.(['im',num2str(i)]) ];
		end        
       
        cal_yet = get(data.handles.cal_text,'String');
        if isfield(data,'cal')
            zoom(dig_fig,'off')
            set(data.handles.pointer_menu,'Value',1)
            set(image_axes,'ButtonDownFcn','mark_points_first_view')
            if isfield(dig_data.temp,'line_1')
                for i = 1:(data.setup.cam_num - 1)
                    line_ax = get(dig_data.temp.(['line_',num2str(i)]),'Parent');
                    line_im = findobj('Type','image','Parent',line_ax);
                    set(line_im,'ButtonDownFcn','mark_points_second_view')
                end
            end
            
        else % If no calibration loaded don't allow to mark points
            zoom(dig_fig,'on')
            set(image_axes,'ButtonDownFcn','')
            set(data.handles.pointer_menu,'Value',2)
        end
        
        % Set appropriate functions/display for coord_plot
        kids = findobj('Parent',dig_data.handles.coord_plot);
        set([dig_data.handles.coord_plot kids'],'ButtonDownFcn','plot_frame_chg')
        
        return
        
    case 'z'
        dig_fig = findobj('Tag','dig_fig');
        dig_data = guidata(dig_fig);

        image_axes = [];
		for i = 1:data.setup.cam_num
            image_axes = [image_axes dig_data.handles.(['im',num2str(i)]) ];
        end
        
        zoom(dig_fig,'on')
        set(image_axes,'ButtonDownFcn','')
        set(data.handles.pointer_menu,'Value',2)
        
        kids = findobj('Parent',dig_data.handles.coord_plot);
        set([dig_data.handles.coord_plot kids'],'ButtonDownFcn','')
        
        return

    case 'u'
        inc_amount = 0.001;
        inc_dec(inc_amount,'x')
        data = guidata(controller); % reload changed data structure
        
    case 'j'
        inc_amount = -0.001;
        inc_dec(inc_amount,'x')
        data = guidata(controller); % reload changed data structure
        
    case 'i'
        inc_amount = 0.001;
        inc_dec(inc_amount,'y')
        data = guidata(controller); % reload changed data structure
        
    case 'k'
        inc_amount = -0.001;
        inc_dec(inc_amount,'y')
        data = guidata(controller); % reload changed data structure
        
    case 'o'
        inc_amount = 0.001;
        inc_dec(inc_amount,'z')
        data = guidata(controller); % reload changed data structure
        
    case 'l'
        inc_amount = -0.001;
        inc_dec(inc_amount,'z')
        data = guidata(controller); % reload changed data structure

    case {'[' '1'}
        frame = str2num(get(data.handles.frame_box,'String'));
        new_frame = frame - 1;
                
        if new_frame < get(data.handles.frame_slider,'Min')
            msgbox('Already at beginning of video sequence.','No more frames','help');
            return
        else
            set(data.handles.frame_box,'String',new_frame)
            update_images
            return
        end
        
    case {']' 'q' '2'}
        frame = str2num(get(data.handles.frame_box,'String'));
        new_frame = frame + 1;
                
        if new_frame > get(data.handles.frame_slider,'Max')
            msgbox('End of video sequence reached.','No more frames','help');
            return
        else
            set(data.handles.frame_box,'String',new_frame)
            update_images
            return
        end
                
        
%     % fine increment of selected axis
%     case 'i'
%         inc_amount = 0.01;
%         inc_dec(inc_amount)
%         data = guidata(controller); % reload changed data structure
%         
%     % coarse decrement of selected axis
%     case 'j'
%         inc_amount = -0.1;
%         inc_dec(inc_amount)
%         data = guidata(controller); % reload changed data structure
%         
%     % fine decrement of selected axis
%     case 'k'
%         inc_amount = -0.01;
%         inc_dec(inc_amount)
%         data = guidata(controller); % reload changed data structure
%         
%     % coarse increment of selected axis
%     case 'l'
%         inc_amount = 0.1;
%         inc_dec(inc_amount)
%         data = guidata(controller); % reload changed data structure
end

guidata(controller, data);

calc_it
plot_it
draw_it


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

    
%     case 'r' % Next Point
%         advance
%         
%     case 'e' % Previous Point
%         advance
%         
%     case '`' % Esc mark point (clear lines)
%         dig_fig = findobj('Tag','dig_fig');
%         dig_data = guidata(dig_fig);
%         
%         % Delete all objects in dig_data.temp
%         points_to_delete = fieldnames(dig_data.temp);
%         for i = 1:length(points_to_delete)
%             delete(dig_data.temp.(points_to_delete{i}))
%             dig_data.temp = rmfield(dig_data.temp,points_to_delete{i});
%         end
%         
%         guidata(dig_fig,dig_data)
%         toggle_pointer_fcn('Mark Points')
%         
%     case 'p' % Clear current point
%         clear_pt('p')
%         
%     case 'o' % Clear current object
%         clear_pt('o')
%         
%     case 'f' % Next frame
%         h = findobj('Tag','frame_box');
%         frame = eval(get(h,'String'));
%         new_frame = frame + 1;
%         
%         set(h,'String',num2str(new_frame))
%         
%         update_images
% 
%     case 'd' % Previous frame
%         h = findobj('Tag','frame_box');
%         frame = eval(get(h,'String'));
%         new_frame = frame - 1;
%         
%         set(h,'String',num2str(new_frame))
%         
%         update_images
%         
%     case 'x' % Set pointer function to Mark Points
%         toggle_pointer_fcn('Mark Points')
%         
%     case 'z' % Set pointer function to Zoom
%         toggle_pointer_fcn('Zoom')
%         
%     case '=' % Save
%         save_data
%         
%     case 'm' % Open object dialog
%         open_object_dialog
%         
%     case 'q' % Reset images (unzoom)
%         reset_button_callback
%         
%     case 'j' % Load images
%         load_images
%         
%     case 'k' % Calibrate
%         load_calibration
%         
%     case 'l' % Load Data
%         load_data
%         
%     case 'v' % Toggle current object visibility
%         h = findobj('Tag','vis_check');
%         cur_val = get(h,'Value');
%         if cur_val == get(h,'Min')
%             new_val = get(h,'Max');
%         else
%             new_val = get(h,'Min');
%         end
%         set(h,'Value',new_val)
%         
%         vis_check_callback
%         
%     case 'b' % Toggle autosave on/off
%         h = findobj('Tag','autosave_check');
%         cur_val = get(h,'Value');
%         if cur_val == get(h,'Min')
%             new_val = get(h,'Max');
%         else
%             new_val = get(h,'Min');
%         end
%         set(h,'Value',new_val)
% 
%         save_auto
