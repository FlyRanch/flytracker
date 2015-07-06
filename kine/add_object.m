% Function ADD_OBJECT
% 
% CALLING FUNCTION: Callback for add_button in obj_fig dialog
% ACTIONS: Sets config_menu to "custom" and opens Edit dialog with "ADD"
%          instead of "okay" button
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 19, 2004 by gwyneth
%
% OBSOLETE!  ADD_OBJECT CALLBACK NOW EDIT_OBJECT_DIALOG

function add_object

% obj_data = guidata(gcf);
% 
% % % Update config_menu to include "custom" and make this the current setting
% % I = strmatch('custom',obj_data.configs,'exact');
% % if isempty(I)
% %     obj_data.configs{end+1} = 'custom';
% % end
% % 
% % set(obj_data.handles.config_menu,'String',obj_data.configs,'Value',length(obj_data.configs))
% 
% % Open Edit Object Dialog
% edit_object_dialog