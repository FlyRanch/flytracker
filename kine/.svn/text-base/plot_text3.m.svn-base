% Function PLOT_TEXT3(ax,x,y,z,st,d,offset_ax)
% 
% CALLING FUNCTION: type_menu_callback(model)
% ACTIONS: Plots text string 5% of z range above or below point entered depending on
%          where in figure point is:
%               ax = axes handle
%               x,y = coordinates of point to label
%               d = distance text offset from point, % of z-range
%               st = string to label point with
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 17, 2004 by gwyneth

function plot_text3(ax,x,y,z,st,d,offset_ax)

% Plot the text below and to the left first, then check to see if goes
% outside figure boundaries

ylim = get(ax,'Ylim');
xlim = get(ax,'XLim');
zlim = get(ax,'ZLim');

ymin = ylim(1);
ymax = ylim(2);
xmin = xlim(1);
xmax = xlim(2);
zmin = zlim(1);
zmax = zlim(2);

% Check to see that original point is in the axes limits
if xmin > x | xmax < x | ymin > y | ymax < y | zmin > z | zmax < z
    return
end

% Plot first go at text
axes(ax)

switch offset_ax
    case '-x'
        text_ht = d/100*(xmax - xmin);
        h = text(x-text_ht, y, z, st);
    case '-y'
        text_ht = d/100*(ymax - ymin);
        h = text(x, y-text_ht, z, st);
    case '-z'
        text_ht = d/100*(zmax - zmin);
        h = text(x, y, z-text_ht, st);
    case '+x'
        text_ht = d/100*(xmax - xmin);
        h = text(x+text_ht, y, z, st);
    case '+y'
        text_ht = d/100*(ymax - ymin);
        h = text(x, y+text_ht, z, st);
    case '+z'
        text_ht = d/100*(zmax - zmin);
        h = text(x, y, z+text_ht, st);
end

% Check to see if extent goes beyond axes limits and accomodate
ext = get(h,'Extent');

if ext(1) + ext(3) > xmax
    set(h,'HorizontalAlignment','right')
end

if ext(2) + ext(4) < ymin
    set(h,'Position',[x y+text_ht])
end
