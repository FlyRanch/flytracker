% Function INITIALIZE_SLIDER
% 
% CALLING FUNCTION: load_images
% ACTIONS: Sets the min, max, and step values of the slider based
%          on how many frames are in the video sequence to load
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 10, 2004 by gwyneth

function initialize_slider

controller = findobj('Tag','controller');
data = guidata(controller);

% Set slider values
set(data.handles.frame_slider,'Min',1);
set(data.handles.frame_slider,'Max',data.images.frames);
set(data.handles.frame_slider,'Value',1);

smallStepLength = 1/(data.images.frames-1);
    if smallStepLength>=1
    smallStepLength=1;
    end
bigStepLength = smallStepLength * 5;
    if bigStepLength>=1
        bigStepLength=smallStepLength*2;
    elseif bigStepLength>=1
        errordlg('something is $#^%@^@$^@$...')
        return      
    end
    
set(data.handles.frame_slider,'SliderStep',[smallStepLength,bigStepLength]);

% Set slider box value to 1
set(data.handles.frame_box,'String',1);
