% Kine picts

function kine_pics(avifilename, frame_range)

avicompress = 'cinepak';
fps = 10;

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');

for frame = frame_range
    
    
    set(data.handles.frame_box,'String',num2str(frame))
    update_images
    
    M(frame) = getframe(dig_fig);
    % pattern_number = frame;
    % n1 = int2str(floor(pattern_number/100));
    % n2 = int2str( floor(mod(pattern_number,100)/10) );
    % n3 = int2str(mod(pattern_number,10));
    % 
    % str = [cur_dir,filesep,avifilename,filesep,avifilename, n1, n2, n3,'.bmp'];
    % 
    % [x,immap] = frame2im(M);
    % imwrite(x, str)
    
end

movie2avi(M,avifilename,'compression',avicompress,'fps', fps)


