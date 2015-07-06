% Function NUM_EDIT_CALLBACK_ARC
% 
% CALLING FUNCTION: Callback for num_edit box in edit_fig when type = arc
% ACTIONS: Displays the appropriate number of pt#_text and pt#_edit
%          uicontrols so that user can enter point names
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function num_edit_callback_arc

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

b_color = get(edit_data.handles.edit_fig,'Color');
num_pts = str2num(get(edit_data.handles.num_edit,'String'));

if num_pts < 3
    set(edit_data.handles.num_edit,'String','4')
    num_pts = 4;
end

arc_text_string = ['To make this arc object, click on the first end point in two frames, then the second end point in two frames, then click on any ',num2str(num_pts),' points alond the arc in the two frames'];
        
% Make arc text and move things around
a = findobj('Tag','arc_text');
if isempty(a)
    ht_below_num = 1;
    text_ht = 10;
    ht_above_buttons = 1;
    
    button_pos = get(edit_data.handles.ok_edit_button,'Position');
    y = button_pos(2) + button_pos(4) + ht_above_buttons;
    y_add = ht_below_num + text_ht + ht_above_buttons;
    
    fig_pos = get(edit_data.handles.edit_fig,'Position');
    fig_pos(2) = fig_pos(2) - y_add; % move down by y_add
    fig_pos(4) = fig_pos(4) + y_add; % make longer by y_add
    set(edit_data.handles.edit_fig,'Position',fig_pos)
    
    for j = 1:length(edit_data.handles.to_move)
        pos = get(edit_data.handles.to_move(j),'Position');
        pos(2) = pos(2) + y_add; % move down by y_add
        set(edit_data.handles.to_move(j),'Position',pos)
    end
    
    arc_text = uicontrol(...
        'Style','text',...
        'Tag','arc_text',...
        'Units','characters',...
        'Position',[10 y 50 text_ht],...
        'BackgroundColor',b_color,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',10,...
        'FontWeight','normal',...
        'HorizontalAlignment','left',...
        'String',arc_text_string...
        );
else
    set(a,'String',arc_text_string)
end

        
