% Function RESET_COORD_SLIDERS
% 
% CALLING FUNCTION: cur_pt_menu_callback, plot_it
% ACTIONS: Resets the range and step value of the coordinate sliders to be
%          +/- the current data range; if no data is entered yet, then
%          defaults to a range of 0-20.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: September 22, 2004 by gwyneth

function reset_coord_sliders

%--------------------------------------------------------------------------
buff = 0.3; % Data range buffer on slider
default_min = 0;
default_max = 1;
%--------------------------------------------------------------------------

controller = findobj('Tag','controller');
data = guidata(controller);

coord = {'x' 'y' 'z'};
cur_obj = get_string(data.handles.cur_object_menu);
cur_pt = get(data.handles.cur_point_menu,'Value');
frame = str2num(get(data.handles.frame_box, 'String'));

for i = 1:3
    
    h_slider = data.handles.([coord{i},'_slider']);
    
    if isfield(data,'kine')
        
        % Find data range
        temp(1,:) = data.kine.(cur_obj).data.coords(cur_pt,i,:);
        [ind_i ind_j v] = find(temp); % at beginning, kine will be empty, v will be []
        
        min_val = min(v);
        max_val = max(v);
        data_range = max_val - min_val;
        
        if isempty(data_range)
            min_val = default_min;
            max_val = default_max;
            data_range = max_val - min_val;
        end
        
        % Slider settings
        if data_range < 10^(-3)
            s_min = min_val - buff*abs(min_val);
            s_max = max_val + buff*abs(max_val);
        else
            s_min = min_val - buff*data_range;
            s_max = max_val + buff*data_range;
        end
        
        % Check for out of range zeros
        cur_val = data.kine.(cur_obj).data.coords(cur_pt,i,frame);
        if cur_val < s_min | cur_val > s_max
            s_min = default_min - buff*default_min;
            s_max = default_max + buff*default_max;
        end
        s_range = s_max - s_min;
        s_step = [ 0.001/s_range 0.01/s_range ];
        set( h_slider, 'Min', s_min, 'Max', s_max, 'SliderStep', s_step )
        
        clear temp v
        
    else
        set(h_slider, 'Min', default_min, 'Max', default_max)
        
    end
    
end
