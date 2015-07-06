% Function MAKE_ALL_VIS
% 
% CALLING FUNCTION: Callback vis_button in obj_fig
% ACTIONS: Sets 'visible' field of all objects to 'on'.  If all 'visible'
%          fields are already set to on, sets them all to 'off'.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function make_all_vis

obj_data = guidata(gcf);

objects = fieldnames(obj_data.kine);

% Check 'visible' field values to see if all 'on'
if isempty(objects) == 0
    for i = 1:length(objects)
        count(i) = strcmp(obj_data.kine.(objects{i}).config.visible,'on');
    end
end

% Set 'visible' fields either all to 'on' or all to 'off'
for i = 1:length(objects)
    if sum(count) == length(objects) % then all visible fields are set to 'on'
        obj_data.kine.(objects{i}).config.visible = 'off';
    else
        obj_data.kine.(objects{i}).config.visible = 'on';
    end
end

guidata(gcf,obj_data)

% Update listbox display
set_listbox_string

% Update config_menu to include "custom" and make this the current setting
add_custom