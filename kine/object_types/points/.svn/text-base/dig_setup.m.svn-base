% Function DIG_SETUP
% 
% CALLING FUNCTION: set_cur_pt_menu (called from obj_ok_callback, etc.)
% ACTIONS: Script for digitizing for type 'POINTS'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 7, 2004 by gwyneth

function dig_setup(cur_obj)

controller = findobj('Tag','controller');
data = guidata(controller);

addpath(data.setup.mfile_path)

% Look through current objects, if don't have 'coords' add zeros
% cur_obj = get_string(data.handles.cur_object_menu);

if isfield(data.kine.(cur_obj),'data') == 0
    data.kine.(cur_obj).data = 'temp';
end

if isfield(data.kine.(cur_obj).data,'coords') == 0
    data.kine.(cur_obj).data.coords = zeros(data.kine.(cur_obj).config.num_pts,3,data.images.frames);
end


guidata(controller,data)

toggle_pointer_fcn('Mark Points')