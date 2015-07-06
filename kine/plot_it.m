% Function PLOT_IT
% 
% CALLING FUNCTION: mark_points_second_view
% ACTIONS: Plots 3D coordinates of the current point in the coordinate plot
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: July 9, 2004 by gwyneth

function plot_it

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Delete old plot lines
%axes(dig_data.handles.coord_plot)
% to_clear = findobj('Tag','coord_plot_line')
% delete(to_clear)
coord_plot_axes = dig_data.handles.coord_plot;
cur_obj_string = get_string(data.handles.cur_object_menu);% Get current point
cur_obj = get_string(data.handles.cur_object_menu);% Get current object
cur_pt_string = get_string(data.handles.cur_point_menu);% Get current point
cur_pt = get(data.handles.cur_point_menu,'Value');% Get current point number
frame = str2num(get(data.handles.frame_box,'String'));

% Only proceed if some coordinates have been entered
if isempty(data.kine.(cur_obj).data.coords) == 0
    % make u-vector (vector of frame numbers)
    u = 1:size(data.kine.(cur_obj).data.coords,3);
        
    if (get(dig_data.handles.angles_check,'Value') == 1) & ...
            (isequal(data.kine.(cur_obj).config.type, 'model')) % box checked and is model, plot angles instead of coordinates

        for i = u
            v(i,1) = rad2deg(data.kine.(cur_obj).data.angles(1,1,i));
            v(i,2) = rad2deg(data.kine.(cur_obj).data.angles(1,2,i));
            v(i,3) = data.kine.(cur_obj).data.params(i);
        end
        title(coord_plot_axes,['Angles of ''',cur_pt_string,''':  phi (RED), theta (GREEN), or alpha (BLUE)'],'FontSize',8,'FontWeight','normal')
    elseif (get(dig_data.handles.angles_check,'Value') == 1) & ...
            (isfield(data.kine.(cur_obj).data, 'eulzyx')) % box checked and is model, plot angles instead of coordinates

        for i = 1:size(data.kine.(cur_obj).data.eulzyx,1)
            %v(i,:) = data.kine.(cur_obj).data.params(i,:);
            v(i,:) = (180/pi)*data.kine.(cur_obj).data.eulzyx(i,:);
        end
        title(coord_plot_axes,['Angles of ''',cur_obj_string,''':  bank (RED), elev (GREEN), or heading (BLUE)'],'FontSize',8,'FontWeight','normal')
    else

        % make v-vector, and plot
        for i = u
            v(i,:) = data.kine.(cur_obj).data.coords(cur_pt,:,i);
        end
        title(coord_plot_axes,['3D coordinates of ''',cur_pt_string,''':  X (RED), Y (GREEN), or Z (BLUE)'],'FontSize',8,'FontWeight','normal')
    end

    if get(dig_data.handles.norm_check,'Value') == 1 % box checked, plot normalized values
        % Get start and stop frames
        start_frame = str2num(get(dig_data.handles.fr_st_edit,'String'));
        stop_frame = str2num(get(dig_data.handles.fr_end_edit,'String'));
	    frames = data.images.frames;                                % get number of frames
        
        if (stop_frame == 0)                                        % if frame range not set then use entire range
            frame_range = 1:frames;
        else
            frame_range = start_frame:stop_frame;
        end
        
        x_min = min(v(frame_range,1));
        y_min = min(v(frame_range,2));
        z_min = min(v(frame_range,3));
        
        v = [v(:,1)-x_min v(:,2)-y_min v(:,3)-z_min];
                
        x_max = max(v(frame_range,1));
        y_max = max(v(frame_range,2));
        z_max = max(v(frame_range,3));
        
        v = [v(:,1)/x_max v(:,2)/y_max v(:,3)/z_max];
        
    end
    size(u);
    size(v);
    plot(coord_plot_axes,u,v,'Tag','coord_plot_line')
            
end

set(dig_data.handles.coord_plot,'NextPlot','add')
ymin = min(get(dig_data.handles.coord_plot,'YLim'));
plot(coord_plot_axes,frame,ymin,'r^','MarkerFaceColor','r','Tag','coord_plot_line')
set(dig_data.handles.coord_plot,'NextPlot','replacechildren')

if get(dig_data.handles.length_check,'Value') == 1 % box checked, show length histogram
    set(dig_data.handles.obj_length,'Visible','on')
    length_ax = dig_data.handles.obj_length;

    col_name = data.kine.(cur_obj).config.color;
    col_rgb = data.colors.(col_name); % look color RGB values up in data.colors

    length_ind = find(data.kine.(cur_obj).data.length); % index of frames with values

    length_mean = mean(data.kine.(cur_obj).data.length(length_ind));
    length_diff = data.kine.(cur_obj).data.length(length_ind) - length_mean;
    cur_length = data.kine.(cur_obj).data.length(frame);
    
    xval = zeros(1,length(length_diff));
    plot(length_ax,xval,length_diff,'.','Color',col_rgb,'Tag','coord_plot_line')
    set(length_ax,'NextPlot','add')
    if cur_length ~= 0 % only plot line if cur_length is not zero
        plot(length_ax,[-1 1],[cur_length cur_length]-length_mean,'r','Tag','coord_plot_line')
    end
    set(length_ax,'NextPlot','replacechildren')
else
    set(dig_data.handles.obj_length,'Visible','off')
end

    
% Update x, y, z slider controls
coord_name = {'x' 'y' 'z'};

for i = 1:3
    
    position = data.kine.(cur_obj).data.coords(cur_pt,i,frame);
    
    h_slider = data.handles.([coord_name{i},'_slider']);
    h_edit = data.handles.([coord_name{i},'_edit']);
    
    % If new point is out of range, reset sliders    
    if position < get(h_slider,'Min') | position > get(h_slider,'Max')
        %set(h_slider,'Min',position-0.15,'Max',position+0.15)
        reset_coord_sliders
    end
    
    set(h_slider,'Value',position)
    set(h_edit,'String',position)
end
% 
% 
%     
% x = data.kine.(cur_object).data.coords(cur_point,1,frame);
% 
% set(data.handles.x_slider,'Value',x)
% set(data.handles.y_slider,'Value',y)
% set(data.handles.z_slider,'Value',z)
% 
% set(data.handles.x_edit,'String',x)
% set(data.handles.y_edit,'String',y)
% set(data.handles.z_edit,'String',z)

set_obj_param_menus