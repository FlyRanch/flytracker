% Function ADVANCE
% 
% CALLING FUNCTION: mark_points_second_view, skip point button
% ACTIONS: Advances current point to the next point (which may be a point
%          in the next object or next frame)
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: March 2, 2004 by gwyneth

function advance

controller = findobj('Tag','controller');
data = guidata(controller);

if gcbo == data.handles.skip_button % the skip button will always advance point, obj, frame
    data.advance = 2;
elseif get(gcf,'CurrentCharacter') == 'r'
    data.advance = 2;
    set(gcf,'CurrentCharacter','y')
elseif get(gcf,'CurrentCharacter') == 'e'
    data.advance = 5;
    set(gcf,'CurrentCharacter','y')
end

cur_point = get(data.handles.cur_point_menu,'Value');

cur_object = get(data.handles.cur_object_menu,'Value');
objects = get(data.handles.cur_object_menu,'String');
cur_object_name = objects{cur_object};

max_points = length(data.kine.(cur_object_name).config.points);
max_objects = length(fieldnames(data.kine));

frame = str2num(get(data.handles.frame_box,'String'));

switch data.advance
    
    case 1 % No advance
        
        return
        
    case 2 % point, object, frame
        
        if cur_point < max_points
            
            new_point = cur_point + 1;
            
            set(data.handles.cur_point_menu,'Value',new_point)
            %draw_it
            cur_pt_menu_callback
            
        elseif cur_point == max_points
            
            if cur_object < max_objects
                
                new_object = cur_object + 1;
                
                set(data.handles.cur_object_menu,'Value',new_object)
                set_cur_pt_menu
                
            elseif cur_object == max_objects
                
                new_frame = frame + 1;
                
                if new_frame > get(data.handles.frame_slider,'Max')
                    ret = questdlg('End of video sequence reached, return to first frame?','No more frames','Yes','No','Yes');
                    switch ret
                        case 'Yes'
                            new_frame = get(data.handles.frame_slider,'Min');
                        case 'No'
                            return
                    end
                end

                set(data.handles.frame_box,'String',new_frame)
                
                set_cur_obj_menu
                update_images
                
            end
        end
        
    case 3 % frame, point, object
        
        new_frame = frame + 1;
        
        if new_frame > get(data.handles.frame_slider,'Max')
            new_frame = get(data.handles.frame_slider,'Min');
            
            if cur_point < max_points
                
                new_point = cur_point + 1;
                set(data.handles.cur_point_menu,'Value',new_point)
                draw_it
                
            elseif cur_point == max_points
                
                if cur_object < max_objects
                    
                    new_object = cur_object + 1;
                    
                    set(data.handles.cur_object_menu,'Value',new_object)
                    set_cur_pt_menu
                    
                elseif cur_object == max_objects
                    
                    ret = questdlg('End of object sequence reached, return to first object?','No more objects','Yes','No','Yes');
                    switch ret
                        case 'Yes'
                            new_object = 1;
                            set(data.handles.cur_object_menu,'Value',new_object)
                            set_cur_pt_menu

                        case 'No'
                            return
                    end

                end
            end
        end
        
        set(data.handles.frame_box,'String',new_frame)
        update_images
        
    case 4 % point, frame, object

        if cur_point < max_points
            
            new_point = cur_point + 1;
            
            set(data.handles.cur_point_menu,'Value',new_point)
            draw_it
            
        elseif cur_point == max_points
            
            new_frame = frame + 1;
            
            if new_frame > get(data.handles.frame_slider,'Max')
                
                new_frame = get(data.handles.frame_slider,'Min');
                
                if cur_object < max_objects
                    
                    new_object = cur_object + 1;
                    
                    set(data.handles.cur_object_menu,'Value',new_object)
                    set_cur_pt_menu
                    
                elseif cur_object == max_objects
                    
                    ret = questdlg('End of object sequence reached, return to first object?','No more objects','Yes','No','Yes');
                    switch ret
                        case 'Yes'
                            new_object = 1;
                            set(data.handles.cur_object_menu,'Value',new_object)
                            set_cur_pt_menu

                        case 'No'
                            return
                    end                    
                end
                
                set(data.handles.frame_box,'String',new_frame)
                
                set_cur_obj_menu
                update_images
                
            end
        end
        
    case 5 % Backwards through point, object, frame
        
        if cur_point > 1
            
            new_point = cur_point - 1;
            
            set(data.handles.cur_point_menu,'Value',new_point)
            %draw_it
            cur_pt_menu_callback
            
        elseif cur_point == 1
            
            if cur_object > 1
                
                new_object = cur_object - 1;
                
                set(data.handles.cur_object_menu,'Value',new_object)
                
                % FROM SET_CUR_PT_MENU
                object = get_string(data.handles.cur_object_menu);
                % Get points for current object
                points = data.kine.(object).config.points;
                if isempty(points)
                    points = {'NONE'};
                end
                
                max_points = length(data.kine.(object).config.points);
                new_point = max_points;
                
                % Set cur_point_menu
                set(data.handles.cur_point_menu,'String',points)
                set(data.handles.cur_point_menu,'Value',new_point)
                
                % Change coordinate plot to be correct point
                calc_it
                plot_it
                draw_it

                
            elseif cur_object == 1
                
                new_frame = frame - 1;
                new_object = max_objects;
                
                if new_frame < get(data.handles.frame_slider,'Min')
                    new_frame = get(data.handles.frame_slider,'Min');
                end
                
                set(data.handles.frame_box,'String',new_frame)
                
                set(data.handles.cur_object_menu,'Value',new_object)
                
                

                % FROM SET_CUR_PT_MENU
                object = get_string(data.handles.cur_object_menu);
                % Get points for current object
                points = data.kine.(object).config.points;
                if isempty(points)
                    points = {'NONE'};
                end
                
                max_points = length(data.kine.(object).config.points);
                new_point = max_points;
                
                % Set cur_point_menu
                set(data.handles.cur_point_menu,'String',points)
                set(data.handles.cur_point_menu,'Value',new_point)
                
                %this happens in update images, so shouldn't need - DEBUG
%                 % Change coordinate plot to be correct point
%                 plot_it
%                 draw_it

                set_obj_param_menus
                
                update_images
                
            end
        end        
        
end % end advance switch
        
        