% Function TYPE_LIST_CALLBACK
% 
% CALLING FUNCTION: Callback for type_list listbox
% ACTIONS: Sets data.kine.(cur_obj).type equal to selected string
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 3, 2004 by gwyneth

function type_list_callback

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu); % Get current object
old_type = data.kine.(cur_obj).config.type;                % Get current type
new_type = get_string(data.handles.type_list);      % Get selected type

if strcmp(old_type,new_type)                        % Do nothing if same type clicked on
    return
end

if (strcmp(old_type,'points') | strcmp(old_type,'segments')) & (strcmp(new_type,'points') | strcmp(new_type,'segments'))       
    proceed = 'Yes';
else
    question = ['Changing the current object to ',new_type,' will delete any data stored for this object.  Do you want to proceed?'];
    title = 'Change object type';
    proceed = questdlg(question, title, 'Yes','No','Yes');    
end

switch proceed
    
    case 'Yes'
        data.kine.(cur_obj).config.type = new_type; 
        % some action will have to be taken here, for now just show in
        % window:
        data.kine.(cur_obj).config.type  %%%DEBUG
    case 'No'
        return
end

guidata(controller,data)