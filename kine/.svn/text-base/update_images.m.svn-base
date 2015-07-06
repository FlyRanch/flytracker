% Function UPDATE_IMAGES
% 
% CALLING FUNCTION: Callback for frame_slider, frame_box; advance
% ACTIONS: Sets the frame box and slider to same values depending on user
%          input and loads the appropriate video frames into all cameras;
%          also runs autosave if data.save.autosave is set to 'on'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 1, 2004 by gwyneth

function update_images

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');

% Find out which object is being used to update images
type = get(gcbo,'Type');
if strcmp(type,'uicontrol')
    obj = get(gcbo,'Style');    % if it was a uicontrol, find out which
else
    obj = 'advance';            % if it's not a uicontrol then it was called by advance (or coord_plot)
end

% Update the other frame-changing object
if strcmp(obj,'slider') == 1
    
    frame = round(get(data.handles.frame_slider,'Value'));
    set(data.handles.frame_box,'String',frame);

else % if edit box set or called by advance.m
    
    % Get value input to frame number box
    frame = eval(get(data.handles.frame_box,'String'));

    % If the entered number is out of the frame range, set box string back to slider value and break
    if frame < get(data.handles.frame_slider,'Min') | frame > get(data.handles.frame_slider,'Max')
        frame = round(get(data.handles.frame_slider,'Value'));
        set(data.handles.frame_box,'String',frame);
        return
    end
   
   % Round to a whole-number index
    frame = round(frame);

    set(data.handles.frame_slider,'Value',frame)
    
end

%display_image

if isfield(data,'cal')  % Only draw on images if calibration has been loaded (data.kine field set up)

    % Set the environment to be correct for the current object type
    object = get_string(data.handles.cur_object_menu);
    type = data.kine.(object).config.type;
    mfile_path = data.setup.mfile_path;
    cd([mfile_path,filesep,'object_types',filesep,type])
    dig_setup(object)
    cd(mfile_path)

    draw_it
    calc_it
    plot_it
end

display_image

% Reset the pointer function to whatever is selected
pointer_choice = get(data.handles.pointer_menu,'Value');
options = get(data.handles.pointer_menu,'String');
toggle_pointer_fcn(options{pointer_choice})

% AUTOSAVE
if strcmp(data.save.autosave,'on')
    save_data
end

figure(dig_fig)