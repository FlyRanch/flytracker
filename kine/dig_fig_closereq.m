% Function DIG_FIG_CLOSEREQ
% 
% CALLING FUNCTION: CloseRequestFcn for dig_fig window
% ACTIONS: deletes dig_fig and deletes handle from data list
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 27, 2004 by gwyneth

function dig_fig_closereq

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');

if dig_fig ~= 0
    question = 'Closing the Digitizing Window will delete all unsaved data.  Do you want to close the window now?';
    title = 'Close Digitizing Window';
    answer = questdlg(question, title, 'Yes','No','Yes');
    
    switch answer
        
        case 'Yes'
            delete(dig_fig)
            
            if isfield(data,'kine') 
                data = rmfield(data,'kine');
            end
            
            if isfield(data,'config')
                data = rmfield(data,'config');
            end
                
            if isfield(data,'images')
                data = rmfield(data,'images');
            end
          
            data.save.pathname = '';
            data.save.filename = '';
            data.save.timestamp = 'Not Saved';
            
            guidata(controller,data)
            update_controller
            
        case 'No'
            return
            
    end
end
