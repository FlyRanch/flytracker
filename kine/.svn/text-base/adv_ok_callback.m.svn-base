% Function ADV_OK_CALLBACK
% 
% CALLING FUNCTION: Callback for OK button in adv_fig
% ACTIONS: Sets data.advance to the value corresponding to the chosen
%          advance order
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: March 4, 2004 by gwyneth

function adv_ok_callback

controller = findobj('Tag','controller');
data = guidata(controller);

advance_fig = findobj('Tag','advance_fig');

% Get handles of all checkboxes
choices = findobj('Parent',advance_fig,'Style','checkbox');

% Find which one is selected
for i = 1:length(choices)
    selected = get(choices(i),'Value');
    
    if selected == get(choices(i),'Max')
        data.advance = get(choices(i),'UserData');
    end
end

guidata(controller,data)
delete(gcf)