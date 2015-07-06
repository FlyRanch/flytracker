% Function UNSPLINE_GUI_CALLBACK
% 
% CALLING FUNCTION: unspline_gui
% ACTIONS: Reads input from unspline_gui and calls unspline_data.m; does
% some error-checking to ensure that the input data values make sense.
% PARENT PROGRAM: Kine_v2_1
% LAST MODIFIED: Jan 08, 2007 by gwyneth

function unspline_gui_callback

boxName = {'start_frame','end_frame','n'};

for i = 1:length(boxName)

    eval([boxName{i},'_handle = findobj(''tag'',''',boxName{i},'_edit'');'])
    eval([boxName{i},' = str2num(get(',boxName{i},'_handle,''String''));'])

    if eval(['isempty(',boxName{i},')'])
        errordlg('Please fill in all values','unspline GUI error')
        return
    end

end

% Must have end_frame > start_frame
% n < (end_frame - start_frame)/2
% no blank values

if start_frame > end_frame
    errordlg('Start Frame cannot be larger than End Frame','unspline GUI error')
    set(start_frame_handle,'String','1')
    return
end

if n > (end_frame - start_frame)/2
    errordlg('N value must be less than half of frame range','unspline GUI error')
    set(n_handle,'String','1')
    return
end

unspline_gui_handle = findobj('tag','unspline_gui');
close(unspline_gui_handle)

unspline_data(start_frame,end_frame,n) % Run unspline data using input values