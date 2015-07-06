% Function SET_TEXT(tag,text,is_field)
% 
% CALLING FUNCTION: load_images
% ACTIONS: Sets the string of the object with the specified tag to the
%          specified text; if the input 'is_field' exists then the 'text'
%          input is interpreted to be the name of a data structure array
%          field.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 10, 2004 by gwyneth


function set_text(tag,text,is_field)

controller = findobj('Tag','controller');
data = guidata(controller);

if exist('is_field')
    controller = findobj('Tag','controller');
    data = guidata(controller);
    
    set(data.handles.(tag),'String',eval(text))
    
else
    
    set(data.handles.(tag),'String',text)
    
end