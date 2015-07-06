% Function SAVE_DATA
% 
% CALLING FUNCTION: Callback for save button, called by save_as
% ACTIONS: Saves variables 'data' and 'timestamp' (created here) into the
%          folder designated by data.setup.save_path
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: March 23, 2004 by gwyneth

function save_data

controller = findobj('Tag','controller');
data = guidata(controller);

while strcmp(data.save.pathname,'')
    data.save.pathname = uigetdir('','Locate the folder in which you would like to save the data');
    if data.save.pathname == 0, return, end
end
while strcmp(data.save.filename,'')
    data.save.filename = inputdlg('Filename: ','Save',1);
    data.save.filename = data.save.filename{1};
    if isempty(data.save.filename), return, end
end
    
if exist(data.save.pathname) == 0
    
    warndlg('Can''t find folder data was previously saved in, try ''Save As...''', 'Save Path Invalid')
    return
end
    
save_name = [data.save.pathname,filesep,data.save.filename];
data.save.timestamp = datestr(now); % current date and time of save

% Save the data
save(save_name, 'data')

% Update the corresponding text
set_text('data_text',['Last save to ''',data.save.filename,''' at ',data.save.timestamp])

guidata(controller,data); % save chosen path and filenames into structure