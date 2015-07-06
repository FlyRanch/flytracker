% Function DELETE_OBJECT
% 
% CALLING FUNCTION: Callback for Delete button edit in obj_fig
% ACTIONS: Deletes object currently selected in listbox
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function delete_object

obj_data = guidata(gcf);

% Get name of currently selected object
obj_name = get_listbox_string;

objects = fieldnames(obj_data.kine);
if isempty(objects)
    warndlg('There are no objects to delete','No object selected')
    return
end

% Query if user is sure the object should be deleted
question = ['Are you sure you want to delete the object:  ',obj_name];
goforit = questdlg(question,'Confirm Delete','Yes','No','Yes');

switch goforit
    
    case 'Yes'
        obj_data.kine = rmfield(obj_data.kine,obj_name);
        guidata(gcf,obj_data)
        
        set_listbox_string
        
        % Update config_menu to include "custom" and make this the current setting
        add_custom
        
    case 'No'
        return
        
end