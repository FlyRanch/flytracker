% export_horizontal_angles.m
% 
% Description:        This function exports angle data (phi, theta, alpha) of all
%                     objects of type 'model' from the data structure to a user
%                     determined file.  The file is in a tab-delimited ASCII format.
%                     
% Arguments:          None.
% 
% Modifies:           None.
% Returns:            None.
% Error handling:     Exits without exporting if user cancells filename selection
%                     dialog box.
%
% Revision history:   07/08/2004  Qing Liu    Initial revision
%                     07/28/2004  Qing Liu    Changed filename; export now draws figure with output plots;
%                                             adjustments to align plots; changed default save filename
%                     08/02/2004  Qing Liu    Changed filename; exports "horizontal world" angles (phi with
%                                             respect to body, theta with respect to real-world horizontal)
%                     

function export_horizontal_angles

[file path no_err] = uiputfile('export_ha.mat', 'Export angles to...');

if (no_err == 0)
    return
end

filename = [path,file];

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

objects = fields(data.kine);                                % get object list
num_objects = length(objects);

frames = data.images.frames;                                % get number of frames

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);                                % retrieve dig_fig guidata

% Get start and stop frames
start_frame = str2num(get(dig_data.handles.fr_st_edit,'String'));
stop_frame = str2num(get(dig_data.handles.fr_end_edit,'String'));

obj_exported = 0;

if (stop_frame == 0)                                        % if frame range not set then use entire range
    frame_range = 1:frames;
else
    frame_range = start_frame:stop_frame;
end

% BEGIN CONSTANT DECLARATIONS
phi = 1; theta = 2;
left_wing=objects{1};body=objects{2};right_wing=objects{3};
head=1;tail=2;
hinge=1;tip=2;
% END CONSTANT DECLARATIONS

coord_diff = data.kine.(body).data.coords(head,:,frame_range) - data.kine.(body).data.coords(tail,:,frame_range);
[theta_b(frame_range) phi_b(frame_range) r_b(frame_range)] = cart2sph(coord_diff(1,1,:), coord_diff(1,2,:), coord_diff(1,3,:));

phi_body = rad2deg(unwrap(theta_b(frame_range))) - 90;     % phi and theta are reversed in MATLAB
theta_body = rad2deg(unwrap(phi_b(frame_range)));
r_body = r_b(frame_range);
    
frames = frame_range(:);
phi_left =      npi2pi(reshape(180 - rad2deg(unwrap(data.kine.(left_wing).data.angles(1,phi,frame_range))), prod(size(frame_range)), 1) + reshape(phi_body, prod(size(frame_range)), 1));
theta_left =    reshape(rad2deg(unwrap(data.kine.(left_wing).data.angles(1,theta,frame_range))),    prod(size(frame_range)), 1);
alpha_left =    reshape(abs(180 - data.kine.(left_wing).data.params(1, frame_range)),            prod(size(frame_range)), 1);
phi_right =     reshape(rad2deg(unwrap(data.kine.(right_wing).data.angles(1,phi,frame_range))), prod(size(frame_range)), 1) - reshape(phi_body, prod(size(frame_range)), 1);
theta_right =   reshape(rad2deg(unwrap(data.kine.(right_wing).data.angles(1,theta,frame_range))),   prod(size(frame_range)), 1);
alpha_right =   reshape(360 - data.kine.(right_wing).data.params(1, frame_range),           prod(size(frame_range)), 1);


save(filename,'frames','phi_left','theta_left','alpha_left','phi_right','theta_right','alpha_right','phi_body','theta_body','r_body')

figure, subplot (3,1,1)
plot (phi_left, 'b')
ylim([-90 90])
hold on
plot (phi_right, 'r')
subplot (3,1,2)
plot (theta_left, 'b')
ylim([-90 90])
hold on
plot (theta_right, 'r')
subplot (3,1,3)
plot (alpha_left, 'b')
ylim([0 180])
hold on
plot (alpha_right, 'r')

figure
subplot (3,1,1)
plot (phi_body)
ylim([-30 30])
subplot (3,1,2)
plot (theta_body)
ylim([0 90])
subplot (3,1,3)
plot (r_body)
ylim([2.75 3.25])


% BEGIN OBSOLETE: the following code exports the angle data to a text file
% data_angles(:, 1) = frame_range(:);                         % output frame numbers to first column
% 
% for obj = 1:num_objects                                     % cycle through objects
%     curr_obj = objects{obj};                                % get object name
%     if (isequal(data.kine.(curr_obj).config.type, 'model')) % if the object is a model, then export angles
%         obj_exported = obj_exported + 1;
%         
%         row_number = 0;
%         for i = frame_range
%             row_number = row_number + 1;
% 
%             data_angles(row_number,(3*(obj_exported-1)+2)) = rad2deg(data.kine.(curr_obj).data.angles(1,2,i));
%             data_angles(row_number,(3*(obj_exported-1)+3)) = rad2deg(data.kine.(curr_obj).data.angles(1,1,i));
%             data_angles(row_number,(3*(obj_exported-1)+4)) = data.kine.(curr_obj).data.params(i);
%         end
%     end
% end
% 
% 
% dlmwrite(filename,data_angles,'\t')
% END OBSOLETE