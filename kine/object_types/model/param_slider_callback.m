% Function PARAM_BUTTON_CALLBACK
% 
% CALLING FUNCTION: callback for +/- buttons in obj control for 'MODEL'
% ACTIONS: change value in param edit box, run param_edit_callback;
%          param_step hard-wired to be 5 (degrees).
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 1, 2004 by gwyneth

function param_button_callback

%--------------------------------------------------------------------------
param_step = 5; %degrees
%--------------------------------------------------------------------------

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);
frame = str2num(get(data.handles.frame_box,'String'));

num = get(gcbo,'UserData');
h_box = findobj('Style','edit','Tag','obj_control','UserData',num);
% old_val = str2num(get(h_box,'String'));

% Add or subtract to get new val, depending on which button pressed
% switch get(gcbo,'String')
%     
%     case '+'
%         new_val = old_val + param_step;
%         
%     case '-'
%         new_val = old_val - param_step;
%         
% end
% 
% % Check limits
% p_min = data.kine.(cur_obj).config.param_array{num,3}(1);
% p_max = data.kine.(cur_obj).config.param_array{num,3}(2);
% 
% if new_val < p_min
%     
%     new_val = p_max - new_val - p_min; % wrap the value around
%     
% elseif new_val > p_max
%     
%     new_val = p_min + new_val - p_max; % wrap the value around
%     
% end

new_val = get(gcbo,'Value');

data.kine.(cur_obj).data.params(num,frame) = new_val;
set(h_box,'String',num2str(round(new_val)))
guidata(controller,data)

obj_program(cur_obj)