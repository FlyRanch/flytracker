% filter_data.m
% 
% DESCRIPTION:      This function filters the x, y, and z values of the
%                   head, tail, and wing hinges. 
%
%                   The current version assumes that these appear as objects 
%                   of the data.kine structure in the following order:
%                   left_wing, body, right_wing, and that within each
%                   object,the hinge appears before the tip and the head 
%                   appears before the tail.
%
%                   The cutoff frequencies are hard-coded.
%             
% ARGUMENTS:        start_frame
%                   end_frame
%             
% LAST MODIFIED:    October 20, 2005 by Doug

function filter_data(start_frame, end_frame)

controller = findobj('Tag','controller');

% Retrieve main guidata
data = guidata(controller);                                 

%set cutoff frequencies
Wn_head = 0.1;
Wn_tail = 0.3;
Wn_wh = 0.3;

% get object list
objects = fields(data.kine);
num_objects = length(objects);

% cycle through objects
for obj = 1:num_objects
    curr_obj = objects{obj};
    num_pts = data.kine.(curr_obj).config.num_pts;
    if obj == 1                 %filter object number one, currently the left wing
        for pt = 1:num_pts
            if pt == 1          %filter the first comonent of the object, currently the wing hinge
                [b,a] = butter(2, Wn_wh)
                for axis = 1:3
                    x = data.kine.(curr_obj).data.coords(pt,axis,:);
                    in_range = find((x >= start_frame) & (x <= end_frame));
                    x = x(in_range);
                    filt_data = filtfilt(b,a,x);
                    for frame = start_frame:end_frame
                        data.kine.(curr_obj).data.coords(pt,axis,frame) = filt_data(frame - start_frame + 1);
                    end
                end
            end
        end
    elseif obj == 2             %filter object number two, currently the body
        for pt = 1:num_pts
            if pt == 1          %filter the first comonent of the object, currently the head
                [b,a] = butter(2, Wn_head)
                for axis = 1:3
                    x = data.kine.(curr_obj).data.coords(pt,axis,:);
                    in_range = find((x >= start_frame) & (x <= end_frame));
                    x = x(in_range);
                    filt_data = filtfilt(b,a,x);
                    for frame = start_frame:end_frame
                        data.kine.(curr_obj).data.coords(pt,axis,frame) = filt_data(frame - start_frame + 1);
                    end
                end
            else                %filter the other comonent of the object, currently the tail
                [b,a] = butter(2, Wn_tail)
                for axis = 1:3
                    x = data.kine.(curr_obj).data.coords(pt,axis,:);
                    in_range = find((x >= start_frame) & (x <= end_frame));
                    x = x(in_range);
                    filt_data = filtfilt(b,a,x);
                    for frame = start_frame:end_frame
                        data.kine.(curr_obj).data.coords(pt,axis,frame) = filt_data(frame - start_frame + 1);
                    end
                end
            end
        end
    elseif obj == 3             %filter object number three, currently the right wing
        for pt = 1:num_pts
            if pt == 1          %filter the first comonent of the object, currently the wing hinge
                [b,a] = butter(2, Wn_wh)
                for axis = 1:3
                    x = data.kine.(curr_obj).data.coords(pt,axis,:);
                    in_range = find((x >= start_frame) & (x <= end_frame));
                    x = x(in_range);
                    filt_data = filtfilt(b,a,x);
                    for frame = start_frame:end_frame
                        data.kine.(curr_obj).data.coords(pt,axis,frame) = filt_data(frame - start_frame + 1);
                    end
                end
            end
        end
    end
end

guidata(controller, data);
