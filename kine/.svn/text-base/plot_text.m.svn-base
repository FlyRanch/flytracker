% Function PLOT_TEXT(ax,x,y,d,st)
% 
% CALLING FUNCTION: type_menu_callback(model)
% ACTIONS: Plots text string 5% of y range above or below point entered depending on
%          where in figure point is:
%               ax = axes handle
%               x,y = coordinates of point to label
%               d = distance text offset from point, % of y-range
%               st = string to label point with
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 17, 2004 by gwyneth

function plot_text(ax,x,y,d,st)

% Plot the text below and to the left first, then check to see if goes
% outside figure boundaries

ylim = get(ax,'Ylim');
xlim = get(ax,'XLim');

ymin = ylim(1);
ymax = ylim(2);
xmin = xlim(1);
xmax = xlim(2);

% Check to see that original point is in the axes limits
if xmin > x | xmax < x | ymin > y | ymax < y
    return
end

% Plot first go at text
axes(ax)
text_ht = d/100*(ymax - ymin);
h = text(x, y-text_ht, st);

% Check to see if extent goes beyond axes limits and accomodate
ext = get(h,'Extent');

if ext(1) + ext(3) > xmax
    set(h,'HorizontalAlignment','right')
end

if ext(2) + ext(4) < ymin
    set(h,'Position',[x y+text_ht])
end
