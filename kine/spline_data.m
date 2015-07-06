% spline_data.m
% 
% Description:    this function performs a spline on the x, y, and z values of points,
%                 as well as recalculates alpha, theta, phi, and radius values of models,
%                 for a user-specified range of frames.
%             
% Arguments:      start_frame
%                 end_frame
%             
% Modifies:       data - structure containing all digitized data
% Returns:        None.
% Error handling: None.

function spline_data(start_frame, end_frame)

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

objects = fields(data.kine);                                % get object list
num_objects = length(objects);

for obj = 1:num_objects                                     % cycle through objects
    curr_obj = objects{obj};                                % get object name
    num_pts = data.kine.(curr_obj).config.num_pts;          % get number of points for that object (g 20040920)
    for pt = 1:num_pts                                      % cycle through points
        for axis = 1:3                                      % cycle through axes
            x = find(data.kine.(curr_obj).data.coords(pt,axis,:) ~= 0);    % find frames with non-zero coordinates
            in_range = find((x >= start_frame) & (x <= end_frame));     % find indices of frames that are within the frame range
            x = x(in_range);                                            % crop to frames that are within range
            if (size(x) < 2)                                % if not enough valid points found, then use start and end frames
                x = [start_frame end_frame];
            end
            y = data.kine.(curr_obj).data.coords(pt,axis,x);            % load values for the frames found above
            xx = start_frame:end_frame;
            yy = spline(x, y, xx);
            for frame = start_frame:end_frame
                data.kine.(curr_obj).data.coords(pt,axis,frame) = yy(frame - start_frame + 1);    % store splined data
            end
            
            % if the object is a model, recalculate theta, phi, and radius
            if (isequal(data.kine.(curr_obj).config.type, 'model')) & (pt == 2)
                for frame = start_frame:end_frame
                    coord_diff = data.kine.(curr_obj).data.coords(2,:,frame) - data.kine.(curr_obj).data.coords(1,:,frame);
                    [theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
                    data.kine.(curr_obj).data.angles(:, :, frame) = [theta phi r];
                end
            end
            
            % if the object is a frame, recalculate theta, phi, and radius
            if (isequal(data.kine.(curr_obj).config.type, 'frame')) & (pt == 2)
                for frame = start_frame:end_frame
                    coord_diff = data.kine.(curr_obj).data.coords(1,:,frame) - data.kine.(curr_obj).data.coords(2,:,frame);
                    [theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
                    data.kine.(curr_obj).data.angles(:, :, frame) = [theta phi r];
                end
            end

            
            if isfield(data.kine.(curr_obj).data,'params')
                
                % only spline alpha once
                if (axis == 1)
                    x = find(data.kine.(curr_obj).data.params(:) ~= 0);    % find frames with non-zero coordinates
                    in_range = find((x >= start_frame) & (x <= end_frame));     % find indices of frames that are within the frame range
                    x = x(in_range);                                            % crop to frames that are within range
                    if (size(x) < 2)                                % if not enough valid points found, then use start and end frames
                        x = [start_frame end_frame];
                    end
                    y = (180/pi)*unwrap( (pi/180)*data.kine.(curr_obj).data.params(x) );
                    xx = start_frame:end_frame;
                    yy = spline(x, y, xx);      
                    
                    wrap_ind = find(yy > 360);
                    yy(wrap_ind) = yy(wrap_ind) - 360;
                    
                    for frame = start_frame:end_frame
                        data.kine.(curr_obj).data.params(frame) = yy(frame - start_frame + 1);    % store splined data
                    end
                end
            end
        end
    end
end

guidata(controller, data);




% % splineData.m
% 
% function splineData
% 
% global KineData;
% head=1;tail=2;rwh=3;rwt=4;lwh=5;lwt=6;rwa=7;lwa=8;
% right=1;left=2; 
% button = KineData.currentButton;
% 
% 
% answer=questdlg('Do you want to spline the data for this button?','Confused scientist dialog','Go ahead','No','no');
% if strcmp(answer,'No')
%    break;
% end
% if ~strcmp(answer,'Go ahead')
%     error('this should not happen...');
% end
% 
% if button <= 6 % all but wing angles
%     % x
%     X=find(KineData.pos.x(button,:) ~= 0);
%     Y=KineData.pos.x(button,X);
%     XX=Y(1):Y(end);
%     XX=X(1):X(end);
%     YY = spline(X,Y,XX);
%     XX=X(1):X(end);
%     KineData.pos.x(button,XX)=YY;
%     
%     % y
%     X=find(KineData.pos.y(button,:) ~= 0);
%     Y=KineData.pos.y(button,X);
%     XX=Y(1):Y(end);
%     XX=X(1):X(end);
%     YY = spline(X,Y,XX);
%     XX=X(1):X(end);
%     KineData.pos.y(button,XX)=YY;
%     
%     % x
%     X=find(KineData.pos.z(button,:) ~= 0);
%     Y=KineData.pos.z(button,X);
%     XX=Y(1):Y(end);
%     XX=X(1):X(end);
%     YY = spline(X,Y,XX);
%     XX=X(1):X(end);
%     KineData.pos.z(button,XX)=YY;
% elseif button == 7
%     X=find(KineData.pos.angles(right,:) ~= 0);
%     Y=KineData.pos.angles(right,X);
%     XX=Y(1):Y(end);
%     XX=X(1):X(end);
%     YY = spline(X,Y,XX);
%     XX=X(1):X(end);
%     KineData.pos.angles(right,XX)=YY;
% elseif button == 8
%     X=find(KineData.pos.angles(left,:) ~= 0);
%     Y=KineData.pos.angles(left,X);
%     XX=Y(1):Y(end);
%     XX=X(1):X(end);
%     YY = spline(X,Y,XX);
%     XX=X(1):X(end);
%     KineData.pos.angles(left,XX)=YY;
%     
% else
%     error('there''s no else here...');
% end
%     
%     updateScreen;
%     