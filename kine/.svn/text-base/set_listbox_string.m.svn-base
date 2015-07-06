% Function SET_LISTBOX_STRING
% 
% CALLING FUNCTION: load_object_config, DeleteFcn for edit_fig
% ACTIONS: looks at obj_data.kine fields and creates string array for
%          object listbox in the object dialog
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 18, 2004 by gwyneth

function set_listbox_string

obj_fig = findobj('Tag','obj_fig');
    if isempty(obj_fig)
        return
    end

obj_data = guidata(obj_fig);

listbox_string = {};

if isfield(obj_data,'kine')
    objects = fieldnames(obj_data.kine);
    if isempty(objects) == 0
        for i = 1:length(objects)
            type = obj_data.kine.(objects{i}).config.type;
            color = obj_data.kine.(objects{i}).config.color;
            vis = obj_data.kine.(objects{i}).config.visible;
            
            if isfield(obj_data.kine.(objects{i}).config,'num_pts')
                num_pts = num2str(obj_data.kine.(objects{i}).config.num_pts);
            else 
                num_pts = '?';
            end
            
            listbox_string{i} = [objects{i},'  (',type,', ',num_pts,' points, ',color,', visible: ',vis,')'];
        end
        
    else
        listbox_string = {'No objects found'};
        
    end
    
else
    listbox_string = {'No objects found'};
    
end

obj_list = findobj('Tag','obj_list');
set(obj_list,'String',listbox_string,'Value',1)