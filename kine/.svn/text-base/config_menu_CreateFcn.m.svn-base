% Function CONFIG_MENU_CREATEFCN
% 
% CALLING FUNCTION: CreateFcn for config_menu in obj_fig; save_config
% ACTIONS: looks at object_configs folder in Kine_v2_0 folder and finds all
%          .mat files, the string of these names is then set as the config_menu
%          string.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 19, 2004 by gwyneth

function config_menu_CreateFcn

controller = findobj('Tag','controller');
data = guidata(controller);

obj_fig = findobj('Tag','obj_fig');
obj_data = guidata(obj_fig);        % Nothing here during creation, but there is for save_config

obj_data.config_path = [data.setup.mfile_path,filesep,'object_configs',filesep];

dir_struct = dir([obj_data.config_path,'*.mat']);     % get names of the .mat files, which are the configs
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

for i = 1:length(sorted_names)
    obj_data.configs{i} = sorted_names{i}(1:end-4);  % get rid of ".mat" at end of file names
end

if isfield(obj_data,'configs') == 0
    obj_data.configs = {'NONE'};
    warndlg('There are no saved configurations to load','No Saved Configurations')
end

config_menu = findobj('Tag','config_menu');
set(config_menu,'String',obj_data.configs)

val = get(config_menu,'Value');            % Keep track of what the config menu is set to
obj_data.cur_config = val;

guidata(obj_fig,obj_data);                % save data structure into figure data structure using GUIDATA

