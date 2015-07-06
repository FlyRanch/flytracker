% Function CONTROLLER_CLOSEREQ
% 
% CALLING FUNCTION: CloseRequestFcn for controller figure
% ACTIONS: Deletes controller and dig_fig windows
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function controller_closereq

controller = findobj('Tag','controller');
dig_fig = findobj('Tag','dig_fig');

if isempty(dig_fig) == 0
    delete(dig_fig)
end

delete(controller)
