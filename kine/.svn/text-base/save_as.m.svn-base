% Function SAVE_AS
% 
% CALLING FUNCTION: Callback for save as menu option
% ACTIONS: Allows user to change file and path name of saved data file
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 1, 2004 by gwyneth

function save_as

controller = findobj('Tag','controller');
data = guidata(controller);

old_path = data.save.pathname;
old_file = data.save.filename;

if strcmp(data.save.pathname,'') & strcmp(data.save.filename,'')
    
    save_data
    
else
    
    data.save.pathname = '';
    data.save.filename = '';
    guidata(controller,data)
    save_data
    
    % if save cancelled, restore old values
    data = guidata(controller);
    
    if strcmp(data.save.pathname,'')
        data.save.pathname = old_path;
    end
    if strcmp(data.save.filename,'')
        data.save.filename = old_file;
    end
    
    guidata(controller,data)
    
end
