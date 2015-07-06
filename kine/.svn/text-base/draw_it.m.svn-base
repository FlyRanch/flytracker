% Function DRAW_IT
% 
% CALLING FUNCTION: mark_points_second_view
% ACTIONS: Plots in all three cameras the coordinates stored in the
%          controller data structure for all objects in the selected frame; Need to
%          add a bit that deals with some objects being made visible to invisible
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: July 6, 2004 by Qing Liu

function draw_it

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Clear all objects with handles in handles structure & handles
% if isfield(dig_data,'objects')%%%DEBUG

if isfield(dig_data,'objects')
    plotted_obj = fieldnames(dig_data.objects);
    
    for i = 1:length(plotted_obj)
        if isfield(dig_data.objects.(plotted_obj{i}),'points')
            
            ind = find(dig_data.objects.(plotted_obj{i}).points);
            non_zero_handles = dig_data.objects.(plotted_obj{i}).points(ind);
            
            h_exist = ishandle(non_zero_handles); % gives vector of 1's and 0's for which ones exist
            ind = find(h_exist);
            non_zero_handles = non_zero_handles(ind);
            
            delete(non_zero_handles)
            dig_data.objects.(plotted_obj{i}) = rmfield(dig_data.objects.(plotted_obj{i}),'points');
        end
    end
    
end

frame = str2num(get(data.handles.frame_box,'String'));

objects = fieldnames(data.kine);
cur_obj_num = get(data.handles.cur_object_menu,'Value');
cur_pt_num = get(data.handles.cur_point_menu,'Value');

for i = 1:length(objects)
    
    cur_object = objects{i};
    col_name = data.kine.(cur_object).config.color;
    col_rgb = data.colors.(col_name); % look color RGB values up in data.colors
    
    if isfield(data.kine.(cur_object),'data')
        num_points = size(data.kine.(cur_object).data.coords,1); %length(data.kine.(cur_object).points);
    else
        num_points = 0;
    end
    
    for j = 1:num_points
        if frame > size(data.kine.(cur_object).data.coords,3)  %%%%DEBUG: change this now that make zero array big enough for all data
            continue
        end
        
        x = data.kine.(cur_object).data.coords(j,1,frame);
        y = data.kine.(cur_object).data.coords(j,2,frame);
        z = data.kine.(cur_object).data.coords(j,3,frame);
        
        
        if x~= 0 & y~=0 %& z~=0
            
            for c = 1:data.setup.cam_num
                
                DLTmat = ['DLT_',num2str(c)];
                if data.setup.cam_num ~= 1
                    [ u(c).pts, v(c).pts ] = dlt_3D_to_2D( data.cal.coeff.(DLTmat), x, y, z );
                else
                    u(c).pts = x;
                    v(c).pts = y;
                end
                
%                axes(dig_data.handles.(['cam',num2str(c)])) 
                cur_ax = dig_data.handles.(['cam',num2str(c)]);
               % cam(c).pt = plot(u(c).pts, v(c).pts, '+', 'MarkerSize',4);
          
               if i == cur_obj_num & j == cur_pt_num
                   mrk = 'o';
                   mrk_size = 3;
               else 
                   mrk = '+';
                   mrk_size = 4;
               end
               
               cam(c).pt = plot(cur_ax, u(c).pts, v(c).pts, mrk,...
                   'MarkerSize', mrk_size, 'Color', col_rgb,...
                   'Visible', data.kine.(cur_object).config.visible);
               
                % Save handles of ploted marks into dig_fig.objects structure
                dig_data.objects.(cur_object).points(j,c) = cam(c).pt;
                
%                 % Set marker color, style appropriately
%                 set(cam(c).pt,'Color',col_rgb,'Visible',data.kine.(cur_object).config.visible)
%                 if i == cur_obj_num & j == cur_pt_num
%                     set(cam(c).pt,'Marker','o','MarkerSize',6)
%                 end
                
            end

            guidata(dig_fig,dig_data)
            dig_data = guidata(dig_fig);
            
        end
        
    end
    
    % Run the obj_program for the correct object type
    type = data.kine.(cur_object).config.type;
    mfile_path = data.setup.mfile_path;
    cd([mfile_path,filesep,'object_types',filesep,type])
    obj_program(cur_object)
    cd(mfile_path)    
    
    % If body-centered views are visible, plot object in bc view
    bc_vis = get(dig_data.handles.body_plot_side,'Visible');
    switch bc_vis
        case 'off'
            continue
        case 'on'
            draw_bcview(cur_object, frame)
    end
    
end
