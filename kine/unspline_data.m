% spline_data.m
% 
% Description:    This function UNSPLINES the data; leaving only every Nth
%                 frame in the specified data range.  Would be nice
%                 if this eventually could be done to individual objects
%                 rather than ALL objects simultaneously.
%             
% Arguments:      start_frame
%                 end_frame
%                 N
%             
% Modifies:       data - structure containing all digitized data
% Returns:        None.
% Error handling: None.
% Modified by:    gwyneth, January 8, 2007

function unspline_data(start_frame, end_frame, N)

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

objects = fields(data.kine);                                % get object list
num_objects = length(objects);

for obj = 1:num_objects                                     % cycle through objects
    curr_obj = objects{obj};                                % get object name
    num_pts = data.kine.(curr_obj).config.num_pts;          % get number of points for that object (g 20040920)

    ind_all = start_frame:end_frame;
    ind_keep = start_frame:N:end_frame;

    if ismember(ind_keep, end_frame) == 0           % add last frame if it's not automatically in range
        ind_keep = [ind_keep end_frame];
    end

    ind_zero = setdiff(ind_all,ind_keep);

    data.kine.(curr_obj).data.coords(:,:,ind_zero) = 0;

    if isfield(data.kine.(curr_obj).data,'params')  % if there are parameters, zero them too
        data.kine.(curr_obj).data.params(:,ind_zero) = 0;
    end

    if isfield(data.kine.(curr_obj).data,'angles')  % if there are angles, zero them too
        data.kine.(curr_obj).data.angles(:,:,ind_zero) = 0;
    end    
end

guidata(controller, data);