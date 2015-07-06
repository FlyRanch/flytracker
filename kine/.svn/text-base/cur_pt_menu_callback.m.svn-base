% Function CUR_PT_MENU_CALLBACK
% 
% CALLING FUNCTION: callback for current point menu
% ACTIONS: plot history of current point, set pointer fuction according to
%          type of point selected
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 3, 2004 by gwyneth


function cur_pt_menu_callback

controller = findobj('Tag','controller');
data = guidata(controller);

%reset_coord_sliders

% plot data stored for this point
calc_it
plot_it
draw_it  % so that current marker is colored red

toggle_pointer_fcn('Mark Points')



% % set the pointer fuction appropriate to the object type and selected point
% cur_obj = get_string(data.handles.cur_object_menu);
% 
% switch data.kine.(cur_obj).config.type
%     
%     case 'points'
%         toggle_pointer_fcn('Mark Points')
%         
%     case 'segments'
%         toggle_pointer_fcn('Mark Points')
%         
%     case 'model'
%         toggle_pointer_fcn('Mark Points')      
%         
%     case 'arc'
%         cur_pt = get_string(data.handles.cur_point_menu);
%         
%         switch cur_pt
%             case 'start'
%                 toggle_pointer_fcn('Mark Points')
%                 
%             case 'end'
%                 toggle_pointer_fcn('Mark Points')
%                 
%             case 'points'
%                 toggle_pointer_fcn('Get Arc Points')
%         end
%         
%     otherwise
%         warndlg('Your object data has not been set correctly','Object type not found')
%         
% end