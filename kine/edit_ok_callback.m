% Function EDIT_OK_CALLBACK
% 
% CALLING FUNCTION: Callback for OK button edit in edit_fig
% ACTIONS: Reads the values of all the uicontrols and outputs these to the
%          obj_fig guidata structure.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function edit_ok_callback

edit_data = guidata(gcf);

controller = findobj('Tag','controller');
data = guidata(controller);

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);


% Get values from all uicontrols using edit_data handles and save into
% obj_data.kine structure
obj_name = get(edit_data.handles.name_edit,'String');
if strcmp(obj_name,'')
    warndlg('You must enter an object name','No object name entered')
    return
end

% Look for spaces in object name (can't have them because have to be field
% names
for i = 1:length(obj_name)
    s = obj_name(i);
    if strcmp(s,' ')
        warndlg('Sorry, can''t have spaces in object names!')
        return
    end
end

% Make sure object has some points
num_pts = str2num(get(edit_data.handles.num_edit,'String'));
if num_pts == 0
    warndlg('Can''t have an object with no points!','No Points')
    return
end

% check for blank edit boxes
edit_boxes = findobj('Parent',obj_fig,'Style','edit');
for i = 1:length(edit_boxes)
    get(edit_boxes(i),'tag')%%DEBUG
    get(edit_boxes(i),'string')
    if strcmp(get(edit_boxes(i),'String'),'');
        warndlg('Please do not leave any boxes blank')
        return
    end
end

% Manipulate fields so objects that already exist are not re-created, just
% modified (if called by the edit button)
objects = fieldnames(obj_data.kine);                                % Get list of objects
if strcmp(edit_data.caller,'edit_button') & length(objects) ~= 0    % Only do this if called by edit button, not add button and if there is an object
    old_obj_name = get_listbox_string;
    if strcmp(obj_name,old_obj_name) == 0                           % Do this if the object name has been changed
        sure = questdlg('You have changed the object name, this will erase any previously saved point data.  Are you sure you want to proceed?');
        if strcmp(sure,'No') | strcmp(sure,'Cancel')                % Warn user that changing object name will clear old data
            return                                                  % Canel ok function if user changes mind
        end
        obj_ind = strmatch(old_obj_name,objects);                   % Find where in the order the object is (by old name)
        obj_data.kine.(obj_name) = obj_data.kine.(old_obj_name);    % Add a new field with new object name equal to old object
        obj_data.kine = rmfield(obj_data.kine,old_obj_name);        % Remove field with old object name
        order_vector = 1:length(objects);                           % Make a vector that will have the order to put the objects in
        order_vector(end) = obj_ind;                                % The newly added field (last field) should have the order of the original object
        after_obj = find(order_vector >= obj_ind);                  % Find all indexes greater than or equal to that of the object
        after_obj = after_obj(1:end-1);                             % Don't include last index which is the object
        for i = after_obj
            order_vector(i) = order_vector(i) + 1;                  % Add one to all of these so none have the same number as the object
        end
        
        new_objects = fieldnames(obj_data.kine);
        for i = 1:length(new_objects)
            new_order(order_vector(i)) = new_objects(i);            % Reorder the fieldnames matrix
        end
        guidata(obj_fig,obj_data)                                   % Save out changes to obj_data
        change_order('obj_fig', 'kine', new_order)                  % Use change order function to rearrange objects to old order
        
    end
end

obj_data = guidata(obj_fig);                                        % Read obj_data in again to get new order

% Can't assign values into subfield of an empty structure, so make
% structure temporarily non-empty then overwrite
if isempty(obj_data.kine)
    obj_data.kine = 'temp';
end

% COLOR
color_options = get(edit_data.handles.color_menu,'String');
choice = get(edit_data.handles.color_menu,'Value');
obj_data.kine.(obj_name).config.color = color_options{choice};

% VISIBLE
vis_check = get(edit_data.handles.vis_box,'Value');
if vis_check == get(edit_data.handles.vis_box,'Max')
    obj_data.kine.(obj_name).config.visible = 'on';
else
    obj_data.kine.(obj_name).config.visible = 'off';
end

% TYPE
type_options = get(edit_data.handles.type_menu,'String');
choice = get(edit_data.handles.type_menu,'Value');
obj_data.kine.(obj_name).config.type = type_options{choice};

% NUM_PTS
obj_data.kine.(obj_name).config.num_pts = str2num(get(edit_data.handles.num_edit,'String'));

% Save out all user-set uicontrol values to the obj_data data structure
guidata(obj_fig,obj_data)   

% POINTS - this depends on type of object, so run script from appropriate
% type folder; assuming that if any fields are left blank this will have
% been caught earlier by edit box checker
mfile_path = data.setup.mfile_path;
type = get_string(edit_data.handles.type_menu);
run([mfile_path,filesep,'object_types',filesep,type,filesep,'edit_ok_callback']);

% Make the current configuration 'custom'
add_custom

% Delete edit object figure
delete(gcf)

