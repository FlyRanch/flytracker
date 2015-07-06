% Function SET_PLOT_VALS
% 
% CALLING FUNCTION: callback for plot edit boxes in dig_fig
% ACTIONS: zooms plot appropriately
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: July 2, 2004 by gwyneth

function set_plot_vals

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Get values to set
xmin = str2num(get(dig_data.handles.fr_st_edit,'String'));
xmax = str2num(get(dig_data.handles.fr_end_edit,'String'));
ymin = str2num(get(dig_data.handles.y_min_edit,'String'));
ymax = str2num(get(dig_data.handles.y_max_edit,'String'));

% Set plot axes with these values
if xmin < xmax
    
    if ymin < ymax
        
        set(dig_data.handles.coord_plot,'XLim',[xmin xmax],'YLim',[ymin ymax])
        
    else
        
        set(dig_data.handles.coord_plot,'XLim',[xmin xmax])
        
        old_val = get(dig_data.handles.coord_plot,'YLim');
        set(dig_data.handles.y_min_edit,'String',num2str(old_val(1)))
        set(dig_data.handles.y_max_edit,'String',num2str(old_val(2)))
        
    end
    
else
    
    if ymin < ymax
        
        set(dig_data.handles.coord_plot,'YLim',[ymin ymax])
        
        old_val = get(dig_data.handles.coord_plot,'XLim');
        set(dig_data.handles.fr_st_edit,'String',num2str(old_val(1)))
        set(dig_data.handles.fr_end_edit,'String',num2str(old_val(2)))
        
    else
        
        old_val = get(dig_data.handles.coord_plot,'YLim');
        set(dig_data.handles.y_min_edit,'String',num2str(old_val(1)))
        set(dig_data.handles.y_max_edit,'String',num2str(old_val(2)))
        
        old_val = get(dig_data.handles.coord_plot,'XLim');
        set(dig_data.handles.fr_st_edit,'String',num2str(old_val(1)))
        set(dig_data.handles.fr_end_edit,'String',num2str(old_val(2)))
        
    end
    
end
