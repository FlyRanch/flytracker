% Function MARK_POINTS_FIRST_VIEW
% 
% CALLING FUNCTION: Button Down function for image axes
% ACTIONS: plot the point on the current axes where the mouse is clicked 
%          and store this point in some kind of data array.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function mark_points_first_view

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);    % get data from dig_fig window

% Do nothing if there is no image file loaded -- don't need to worry about
% this because the callback is from the image, not the axes, so no callback
% should be activated if no image is loaded.

% Retrieve relevant info:
coords = get(gca,'CurrentPoint');    
x = coords(1,1);
y = coords(1,2);

cam = get(gca,'Tag');
cam = str2num(cam(end));

available_cams = [1:data.setup.cam_num];
others = find([available_cams] ~= cam);

% Plot marked point
if isfield(dig_data.temp,'point')
    set(dig_data.temp.point,'XData',x,'YData',y,'ButtonDownFcn','mark_points_first_view')
else
    dig_data.temp.point = plot(x,y,'r+','Tag',num2str(cam)); % The tag of the point is the number of the cam used
    guidata(dig_fig,dig_data);
    toggle_pointer_fcn('find_3d',others)
end

plot_lines(x,y,cam,others)

% Set Button down function so that can click over lines for determining second point
dig_data = guidata(dig_fig);

for i = others
    set(dig_data.temp.(['line_',num2str(i)]),'ButtonDownFcn','mark_points_second_view')
end


if isempty(others) % Case where have only one camera -- just save 2D points
    
    H = [ x y 0];
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
    
end