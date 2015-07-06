% Function MARK_POINTS_SECOND_VIEW
% 
% CALLING FUNCTION: ButtonDownFcn for 'other' cameras after point in first
%                   view selected during digitizing
% ACTIONS: The Button Down callback for the 'other' images (the ones where 
%          no point selected); uses info from first click and current
%          click to determine 3D coordinates.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function mark_points_second_view

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Retrieve relevant info:
    % First point coordinates
    x_cam_1 = get(dig_data.temp.point,'XData');
    y_cam_1 = get(dig_data.temp.point,'YData');
    
    % DLT coeffs corresponding to camera first point is from
    cam_1 = get(dig_data.temp.point,'Tag');
    DLT_cam_1 = data.cal.coeff.(['DLT_',cam_1]);
        
    % Second point coordinates
    coords = get(gca,'CurrentPoint');    
    x_cam_2 = coords(1,1);
    y_cam_2 = coords(1,2);
    
    % DLT coeffs corresponding to camera second point is from
    cam_2 = get(gca,'Tag');
    cam_2 = cam_2(end);
    DLT_cam_2 = data.cal.coeff.(['DLT_',cam_2]);

% Reconfu with first and second points clicked on
A = [DLT_cam_1, DLT_cam_2];
L = [x_cam_1, y_cam_1, x_cam_2, y_cam_2];

[H] = reconfu(A,L);

% Save coordinates in controller data structure
store(H)

% Delete all objects in dig_data.temp
points_to_delete = fieldnames(dig_data.temp);
for i = 1:length(points_to_delete)
    delete(dig_data.temp.(points_to_delete{i}))
    dig_data.temp = rmfield(dig_data.temp,points_to_delete{i});
end

guidata(dig_fig,dig_data)

draw_it
advance
calc_it
plot_it

toggle_pointer_fcn('Mark Points')
    
% -------------------------------------------------------------------------
% Note, the calling syntax for reconfu is:

%function [H] = reconfu(A,L)
% Description:  Reconstruction of 3D coordinates with the use local (camera
%               coordinates and the DLT coefficients for the n cameras).
% Input:        - A  file containing DLT coefficients of the n cameras
%                    [a1cam1,a1cam2...;a2cam1...]
%               - L  camera coordinates of points
%                    [xcam1,ycam1,xcam2,ycam2...;same at time 2]
% Output:       - H  global coordinates, residuals, cameras used
%                    [Xt1,Yt1,Zt1,residt1,cams_used@t1...; same for t2]
