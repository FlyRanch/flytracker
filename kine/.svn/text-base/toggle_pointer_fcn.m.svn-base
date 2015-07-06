% Function TOGGLE_POINTER_FCN(fcn, varargin)
% 
% CALLING FUNCTION: Callback for pointer_menu; load_images, obj_ok_callback, draw_it
% ACTIONS: Toggles between pointer functions, options:
%               fcn = 'Mark Points'
%                     'Zoom'
%                     'Clear Points'
%                     'Clear Object'
%                     'find_3d'
%          pointer task is initiated when the pointer is  
%          clicked in any of the image axes (cam1, cam2, cam3).
%          Inputs after the first are for use with 'find_3d', and specify
%          which camera images should be set to mark_points_second_view
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function toggle_pointer_fcn(fcn,varargin)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

if strcmp(fcn,'menu_in')
    choice = get(data.handles.pointer_menu,'Value');
    options = get(data.handles.pointer_menu,'String');
    fcn = options{choice};
end
    
% Get the handles of all the images
image_axes = [];
for i = 1:data.setup.cam_num
    image_axes = [image_axes dig_data.handles.(['im',num2str(i)]) ];
end

% Get the handles of all the drawn objects and points
all_lines = findobj('type','line');
plot_lines = findobj('tag','coord_plot_line');
obj_lines = setdiff(all_lines, plot_lines);

switch fcn
    
    case 'Mark Points'
        
        cal_yet = get(data.handles.cal_text,'String');
        if isfield(data,'cal')
            zoom(dig_fig,'off')
            set(data.handles.pointer_menu,'Value',1)
            set(image_axes,'ButtonDownFcn','mark_points_first_view')
%             if isfield(dig_data.temp,'line_1')% should be redundant with next command
%                 for i = 1:(data.setup.cam_num - 1)
%                     line_ax = get(dig_data.temp.(['line_',num2str(i)]),'Parent');
%                     line_im = findobj('Type','image','Parent',line_ax);
%                     set(line_im,'ButtonDownFcn','mark_points_second_view')
%                 end
%             end
            
            set(obj_lines,'ButtonDownFcn','mark_points_first_view')
            
        else % If no calibration loaded don't allow to mark points
            zoom(dig_fig,'on')
            set(image_axes,'ButtonDownFcn','')
            set(obj_lines,'ButtonDownFcn','')
            set(data.handles.pointer_menu,'Value',2)
        end
        
        % Set appropriate functions/display for coord_plot
        kids = findobj('Parent',dig_data.handles.coord_plot);
        set([dig_data.handles.coord_plot kids'],'ButtonDownFcn','plot_frame_chg')
        

    case 'Zoom'
        
        zoom(dig_fig,'on')
        set(image_axes,'ButtonDownFcn','')
        set(obj_lines,'ButtonDownFcn','')
        set(data.handles.pointer_menu,'Value',2)
        
        kids = findobj('Parent',dig_data.handles.coord_plot);
        set([dig_data.handles.coord_plot kids'],'ButtonDownFcn','')
        set(dig_data.handles.dig_fig,'Pointer','arrow');

        
    case 'Clear Point'
        
        zoom off
        set(image_axes,'ButtonDownFcn','clear_pt')
        set(obj_lines,'ButtonDownFcn','')
        
    case 'Clear Object'
        
        zoom off
        set(image_axes,'ButtonDownFcn','clear_pt')
        set(obj_lines,'ButtonDownFcn','')
        
    case 'find_3d'
        
        zoom off
        
        % Should use the first variable input after 'fcn'
        % Correct input for varargin would be a vector with cam numbers
        % that should have button down set to 'second_point'
%         set(image_axes(varargin{1}),'ButtonDownFcn','mark_points_second_view') %'other' images
        
        % get correct obj_lines handles
        h_inc = [];
        for i = 1:length(varargin{1})
            
            h_ax = findobj('tag',['cam',num2str(varargin{1}(i))]);
%             h_im = findobj('tag',['im',num2str(varargin{1})]);
            h_all = get(h_ax,'children');
            
            h_inc = [h_inc; h_all]; %setdiff(h_all,h_im)]
        end
                        
        set(h_inc,'ButtonDownFcn','mark_points_second_view')
        
        
end
