% Function ADV_CHECKBOX_CALLBACK
% 
% CALLING FUNCTION: Callback for Advance checkboxes
% ACTIONS: Allow only one checkbox to be selected at a time
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: March 4, 2004 by gwyneth

function adv_checkbox_callback

advance_fig = findobj('Tag','advance_fig');

% Get handles of all option checkboxes
choices = findobj('Parent',advance_fig,'Style','checkbox');

% Set all the ones not called to min value

ind = find(choices~=gcbo);
not_selected = choices(ind);
set(not_selected,'Value',get(not_selected(1),'Min'))