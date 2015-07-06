% Function EDIT_OK_CALLBACK
% 
% CALLING FUNCTION: Kine_v2_0
% ACTIONS: Edit dialog okay button callback for type 'MODEL'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 17, 2004 by gwyneth

function edit_ok_callback

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

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

% add model fields to obj_data
obj_data.kine.(obj).config.points = anchor_array(:,1);
obj_data.kine.(obj).config.model_name = model_name;
obj_data.kine.(obj).config.model_coords = coords;
obj_data.kine.(obj).config.anchor_array = anchor_array;
% obj_data.kine.(obj).config.param_array = param_array; NO PARAMETERS YET

guidata(obj_fig,obj_data)
guidata(controller, data)