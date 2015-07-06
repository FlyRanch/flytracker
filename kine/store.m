% Function STORE(coords)
% 
% CALLING FUNCTION: mark_points_second_view
% ACTIONS: Stores current point information to correct part of  
%          controller data structure array ; takes "coords" from 
%          reconfu output (so takes only the first 3 values of 
%          coords)
% PARENT PROGRAM: Kine_v2_0
% REVISION HISTORY: March 1, 2004 by gwyneth
%                   June 19, 2004 by Qing Liu - added angle calculation
%                       code for objects of type 'model'

function store(coords)

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

coords = coords(1:3);                                       % Take first three values of coords matrix because 
                                                            % these are the x,y,z, values from reconfu
                                                            
% Get current object,  point, and frame from menus
choice = get(data.handles.cur_object_menu,'Value');
objects = get(data.handles.cur_object_menu,'String');
cur_object = objects{choice};                               % Note that while the cur_object is retrieved as a string
                                                            % because it will be used to retrieve a structure field, the
cur_point = get(data.handles.cur_point_menu,'Value');       % cur_point is retrieved as a number because it will be used
                                                            % to index the rows of a matrix - this number should be okay without
frame = str2num(get(data.handles.frame_box,'String'));      % verification because the menu is made by retrieving the points
                                                            % string, so the name of the point will be data.kine.(object).points{cur_point}
                                                            
% Save 3D coords into appropriate data structure
data.kine.(cur_object).data.coords(cur_point,:,frame) = coords;  % This should work even if the field (cur_object) does not exist

% If point is one in a model then calculate the angles associated with the
% two points and save
if isequal(data.kine.(cur_object).config.type, 'model') | isequal(data.kine.(cur_object).config.type, 'frame') 
        coord_diff = data.kine.(cur_object).data.coords(2,:,frame) - data.kine.(cur_object).data.coords(1,:,frame);
        [theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
        data.kine.(cur_object).data.angles(:, :, frame) = [theta phi r];
end

% Update body centered coordinates and angles
% calc_angles_stroke_plane_centered(frame)  COMMENTED OUT BY QING, PUT BACK
% EVENTUALLY

guidata(controller,data);