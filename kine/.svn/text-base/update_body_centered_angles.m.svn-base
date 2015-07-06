function update_body_centered_angles

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

frames = data.images.frames;                                % get number of frames

for frame = 1:frames
    calc_angles_stroke_plane_centered(frame);
end

