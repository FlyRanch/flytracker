% Function EDIT_OK_CALLBACK
% 
% CALLING FUNCTION: Kine_v2_0
% ACTIONS: Edit dialog okay button callback for type 'FRAME'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: October 2, 2004 by gwyneth

function edit_ok_callback

controller = findobj('Tag','controller');
data = guidata(controller);

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

if isfield(edit_data,'model')
    model_name = edit_data.model.name;
    model_dir = edit_data.model.dir;
else
    model_name = edit_data.model_name;
    model_dir = [data.setup.mfile_path,filesep,'models'];
end


obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);

% get object name
obj = get(edit_data.handles.name_edit,'String');

% get model name
h = findobj('Tag','model name text');
model_name = get(h,'String');
model_name = model_name(17:end);

% load model
cd(model_dir)
load(model_name)
cd(data.setup.mfile_path)

% Get info from obj/point menus to make match array
for i = 1:length(match_points)

    h_ob = findobj('Tag',['ob_menu_',num2str(i)]);
    h_pt = findobj('Tag',['pt_menu_',num2str(i)]);
    
    match_array(i,:) = { match_points(i) get_string(h_ob) get_string(h_pt) };
    
end

% add model fields to obj_data
obj_data.kine.(obj).config.points = anchor_array(:,1);
obj_data.kine.(obj).config.model_name = model_name;
obj_data.kine.(obj).config.model_coords = coords;
obj_data.kine.(obj).config.anchor_array = anchor_array;
obj_data.kine.(obj).config.param_array = param_array;
obj_data.kine.(obj).config.match_points = match_points;
obj_data.kine.(obj).config.match_array = match_array;

guidata(obj_fig,obj_data)
guidata(controller, data)