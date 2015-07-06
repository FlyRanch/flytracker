% Function LOAD_IMAGES
% 
% CALLING FUNCTION: load_button callback
% ACTIONS: Opens Digitizing window, initializes slider, and loads images for first frame
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 3, 2004 by gwyneth

function load_images

dig_fig_closereq                % Close any open dig fig windows

proceed = findobj('Tag','dig_fig');
if proceed ~= 0
    return                      % Only proceed if the user agreed to close the current dig_fig window
end

success = get_image_file_info;  % Query user for image file location/name
if success == 1
    initialize_slider           % Set frame slider parameters based on number of frames detected
    open_dig_fig                % Creates Digitizing window (contains dig_fig_closereq, but dig_fig should be closed, so shouldn't run)
    display_image               % Display the images in the three camera axes
    update_controller           % update controller display
end