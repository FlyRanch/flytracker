% Function CROSSHAIR_HOTSPOTS
% 
% CALLING FUNCTION: WindowButtonMotionFcn for dig_fig (see open_dig_fig.m)
% ACTIONS: Changes pointer to a crosshair whenever it moves over certain "hotspots" 
%          in the figure window.  For current purposes these hotspots are the image 
%          axes with handles:
%                           data.handles.cam1
%                           data.handles.cam2
%                           data.handles.cam3
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 30, 2004 by gwyneth

function crosshair_hotspots

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Get current pointer location in dig_fig window
cur_point = get(dig_fig,'CurrentPoint');

%--------------------------------------------------------------------------
%CROSSHAIR HOTSPOTS

% Get handles of all axes hotspots
for i = 1:data.setup.cam_num
    hotspot = get(dig_data.handles.(['cam',num2str(i)]),'Position');
    x_range(i,:) = [hotspot(1) hotspot(1)+hotspot(3)];
    y_range(i,:) = [hotspot(2) hotspot(2)+hotspot(4)];
end

% Determine if pointer is in any of the hotspots
for i = 1:data.setup.cam_num
    if x_range(i,1) <= cur_point(1) & x_range(i,2) >= cur_point(1)
        x_match(i) = 1;
    else
        x_match(i) = 0;
    end
        if y_range(i,1) <= cur_point(2) & y_range(i,2) >= cur_point(2)
        y_match(i) = 1;
    else
        y_match(i) = 0;
    end
end

if dot(x_match,y_match)
    set(dig_data.handles.dig_fig,'Pointer','crosshair');
else
    set(dig_data.handles.dig_fig,'Pointer','arrow');
end

%--------------------------------------------------------------------------
%FRAME SELECT HOTSPOTS

% Get coordinates of the plot axes
pos = get(dig_data.handles.coord_plot,'Position');

% Set hotspot around the bottom axis
del_h = 0.005;

w_range = [pos(1) pos(1)+pos(3)];
h_range = [pos(2)-del_h pos(2)+pos(4)];

if w_range(1) <= cur_point(1) & w_range(2) >= cur_point(1)...
        & h_range(1) <= cur_point(2) & h_range(2) >= cur_point(2)
    set(dig_data.handles.dig_fig,'Pointer','circle');
elseif strmatch(get(dig_data.handles.dig_fig,'Pointer'),'crosshair') == 0
    set(dig_data.handles.dig_fig,'Pointer','arrow');
end
