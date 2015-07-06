% Function PLOT_LINES(x,y,cam,others)
% 
% CALLING FUNCTION: Button Down function for image axes
% ACTIONS: plot lines in other camera views that correspond to
%          the point clicked on in the current camera view; inputs are:
%               - x,y     2D coordinates of point in view 1
%               - cam     number of cam corresponding to view 1
%               - others  numbers of other cams where lines are to be drawn
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function plot_lines(x,y,cam,others)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Assign proper labels
point_cam = dig_data.handles.(['cam',num2str(cam)]);
point_cam_DLT = data.cal.coeff.(['DLT_',num2str(cam)]);

for i = others
    
    % Assign proper labels
    line(i).cam = dig_data.handles.(['cam',num2str(i)]); 
    line(i).im = dig_data.handles.(['im',num2str(i)]); 
    line(i).DLT = data.cal.coeff.(['DLT_',num2str(i)]);
    
    % Detect image resolution for current "other" camera
    xres = get(line(i).im,'XData');
    yres = get(line(i).im,'YData');
    
    % Get line coordinates
    [u(i).pts, v(i).pts] = im_pt_2_im_line( x, y, point_cam_DLT, line(i).DLT, [xres yres], 10 );
    axes(line(i).cam)
    
    line(i).name = ['line_',num2str(i)];
    if isfield(dig_data.temp,line(i).name)
        set(dig_data.temp.(line(i).name),'XData',u(i).pts,'YData',v(i).pts)
    else
        dig_data.temp.(line(i).name) = plot(u(i).pts, v(i).pts,'Color','r','Tag',num2str(i));
    end
    
end

guidata(dig_fig,dig_data);