% Function SAVE_AUTO
% 
% CALLING FUNCTION: Callback for autosave menu option
% ACTIONS: Toggles autosave setting; when on, autosave will automatically
%          run save_data.m whenever the frame is changed (in update_images)
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 1, 2004 by gwyneth

function save_auto

controller = findobj('Tag','controller');
data = guidata(controller);

state = get(data.handles.autosave_check,'Value');
check = get(data.handles.autosave_check,'Max');
no_check = get(data.handles.autosave_check,'Min');

switch state
    
    case check 
        data.save.autosave = 'on';
        
    case no_check % if it's off, turn it on
        data.save.autosave = 'off';
        
end

guidata(controller,data)
