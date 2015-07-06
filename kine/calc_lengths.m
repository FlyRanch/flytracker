% Function CALC_LENGTHS
% 
% CALLING FUNCTION: calculate_it
% ACTIONS: Calculates lengths of all objects in the frame that:
%            (1) have been digitized (no zero elements)
%            (2) have at least two points associated with them
%          Saves the length to data.kine.(obj).data.length.
% PARENT PROGRAM: Kine_v2_1
% LAST MODIFIED: November 20, 2006 by gwyneth

function calc_lengths

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

if get(dig_data.handles.length_check,'Value') == 1 % box checked, show length histogram
    frame = str2num(get(data.handles.frame_box,'String'));
    objects = fieldnames(data.kine);

    for i = 1:length(objects)

        cur_object = objects{i};
        col_name = data.kine.(cur_object).config.color;
        col_rgb = data.colors.(col_name); % look color RGB values up in data.colors
        
        if data.kine.(cur_object).data.coords(1,1,frame) == 0 | data.kine.(cur_object).data.coords(2,1,frame) == 0
            L = 0;
        else

            for j = 1:2 % BASE LENGTH ON DISTANCE BETWEEN FIRST TWO POINTS
                if frame > size(data.kine.(cur_object).data.coords,3)  %%%%DEBUG: change this now that make zero array big enough for all data
                    continue
                end


                x(j) = data.kine.(cur_object).data.coords(j,1,frame);
                y(j) = data.kine.(cur_object).data.coords(j,2,frame);
                z(j) = data.kine.(cur_object).data.coords(j,3,frame);
            end

            L = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 + (z(1)-z(2))^2 );
        end

        data.kine.(cur_object).data.length(frame) = L;

    end

    guidata(controller,data)
end
