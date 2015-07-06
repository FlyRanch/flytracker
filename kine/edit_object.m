% Function EDIT_OBJECT
% 
% CALLING FUNCTION: Callback for Edit button edit in obj_fig
% ACTIONS: Opens edit object dialog then sets values of uicontrols
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function edit_object

controller = findobj('Tag','controller');
data = guidata(controller);

% Create Edit Object Dialog
edit_object_dialog

% Set value of uicontrols based on currently selected object
obj_fig = findobj('Tag','obj_fig');
edit_fig = findobj('Tag','edit_fig');

obj_data = guidata(obj_fig);
edit_data = guidata(edit_fig);

% If there are no objects in this configuration, return without setting
% uicontrol values
if isempty(obj_data.kine)
    return
end

% Get selected object name from listbox
obj_name = get_listbox_string;

% Set object name in edit dialog figure
set(edit_data.handles.name_edit,'String',obj_name)

% Set object color
color_options = get(edit_data.handles.color_menu,'String');
cur_color = obj_data.kine.(obj_name).config.color;
I = strmatch(cur_color,color_options,'exact');
set(edit_data.handles.color_menu,'Value',I)

% Set visibility
on = get(edit_data.handles.vis_box,'Max');
off = get(edit_data.handles.vis_box,'Min');

if strcmp(obj_data.kine.(obj_name).config.visible,'on')
    set(edit_data.handles.vis_box,'Value',on)
else
    set(edit_data.handles.vis_box,'Value',off)
end

% Set object type
type_options = get(edit_data.handles.type_menu,'String');
cur_type = obj_data.kine.(obj_name).config.type;
I = strmatch(cur_type,type_options,'exact');
set(edit_data.handles.type_menu,'Value',I)

% Set number of points
num_pts = obj_data.kine.(obj_name).config.num_pts;
set(edit_data.handles.num_edit,'String',num2str(num_pts))

% Set point field according to type
mfile_path = data.setup.mfile_path;
run([mfile_path,filesep,'object_types',filesep,cur_type,filesep,'edit_object'])

