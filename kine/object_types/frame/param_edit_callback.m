% Function PARAM_EDIT_CALLBACK
% 
% CALLING FUNCTION: callback for edit boxes in obj control for 'MODEL'
% ACTIONS: Enter current edit box value into appropriate parameter array
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 1, 2004 by gwyneth

function param_edit_callback

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);
frame = str2num(get(data.handles.frame_box,'String'));
num = get(gcbo,'UserData');
val = str2num(get(gcbo,'String'));
h_slider = findobj('Tag','obj_control','Style','slider','UserData',num);
%get(h_slider)

% Check limits
p_min = data.kine.(cur_obj).config.param_array{num,3}(1);
p_max = data.kine.(cur_obj).config.param_array{num,3}(2);

if val >= p_min & val <= p_max
    
    data.kine.(cur_obj).data.params(num,frame) = val;
    set(h_slider,'Value',val)
    
else
    
    warndlg(['The value for parameter ',data.kine.(cur_obj).config.param_array{num,1},' must be between ',num2str(p_min),' and ',num2str(p_max)])
    set(gcbo,'String',data.kine.(cur_obj).data.params(num,frame))
    
end

guidata(controller,data)
obj_program(cur_obj)
