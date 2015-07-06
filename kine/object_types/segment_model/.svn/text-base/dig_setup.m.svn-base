% Function DIG_SETUP
% 
% CALLING FUNCTION: Kine_v2_0
% ACTIONS: Script for digitizing for type 'MODEL'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 7, 2004 by gwyneth

function dig_setup(cur_obj)

controller = findobj('Tag','controller');
data = guidata(controller);

addpath(data.setup.mfile_path)

% Look at current object, if don't have 'coords' and 'params' add zeros
% cur_obj = get_string(data.handles.cur_object_menu);

if isfield(data.kine.(cur_obj),'data') == 0
    data.kine.(cur_obj).data = 'temp';
end

if isfield(data.kine.(cur_obj).data,'coords') == 0
    data.kine.(cur_obj).data.coords = zeros(data.kine.(cur_obj).config.num_pts,3,data.images.frames);
end

% NO PARAMETERS YET
% if isfield(data.kine.(cur_obj).data,'params') == 0
%     data.kine.(cur_obj).data.params = zeros(size(data.kine.(cur_obj).config.param_array,1),data.images.frames);
% end
% 
% if isfield(data.kine.(cur_obj).data,'angles') == 0
%     % Make data.kine and data.config
%     data.kine.(cur_obj).data.angles = zeros(1,3,data.images.frames);   % if the object is a model, initialize the angles array
% end

guidata(controller,data)
toggle_pointer_fcn('Mark Points')

% Delete any existing object controls
to_delete = findobj('Tag','obj_control');
delete(to_delete)

% % Check to see if all the anchor points of the current object are digitized
% % for the current frame
% % cur_obj = get_string(data.handles.cur_object_menu);
% frame = str2num(get(data.handles.frame_box,'String'));
% 
% undone = find(data.kine.(cur_obj).data.coords(:,:,frame) == 0);
% if isempty(undone) == 0
% 
%     figure(controller)
% 
%     uicontrol(...
%         'Style','pushbutton',...
%         'Tag','obj_control',...
%         'Units','normalized',...
%         'Position',[0.2 0.05 0.6 0.1],...
%         'BackgroundColor',data.colors.button_color,...
%         'FontName','Helvetica',...
%         'FontUnits','points',...
%         'FontSize',11,...
%         'FontWeight','bold',...
%         'String','FIT MODEL',...
%         'ToolTipString','Click here to fit the model and set model parameters',...
%         'Visible','on',...
%         'Callback','model_params'...
%         );
%     
% else
% 
%     model_params
%     
% end

obj_program(cur_obj)