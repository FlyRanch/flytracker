% Function NUM_EDIT_CALLBACK
% 
% CALLING FUNCTION: Callback for num_edit edit box in edit_fig
% ACTIONS: Displays the appropriate number of pt#_text and pt#_edit
%          uicontrols so that user can enter point names
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 20, 2004 by gwyneth

function num_edit_callback

%--------------------------------------------------------------------------
edit_ht = 1.75;
text_ht = edit_ht - 0.25;
space_ht = 1;
base_ht = 5;

x_text = 23;
x_edit = x_text + 14;
%--------------------------------------------------------------------------

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

b_color = get(edit_data.handles.edit_fig,'Color');
num_pts = str2num(get(edit_data.handles.num_edit,'String'));

type_options = get(edit_data.handles.type_menu,'String');
choice = get(edit_data.handles.type_menu,'Value');
type = type_options{choice};


% If num_pts entered is less than 1 or greater than 10, set to 1 or 100
if num_pts == 0
    num_pts = 1;
    set(edit_data.handles.num_edit,'String',num2str(num_pts))
elseif num_pts > 100
    num_pts = 100;
    set(edit_data.handles.num_edit,'String',num2str(num_pts))
end

% clear any existing point texts and edits
to_clear = findobj('Tag','pt_text');
delete(to_clear)

to_clear = findobj('Tag','pt_edit');
delete(to_clear)

% make appropriate text and edit boxes, resize dialog window
old_n = length(to_clear);
old_y_add = old_n*(edit_ht + space_ht); % correct for previous movement

y_add = num_pts*(edit_ht + space_ht);   % resize according to number of points

fig_pos = get(edit_data.handles.edit_fig,'Position');
fig_pos(2) = fig_pos(2) + old_y_add - y_add; % move down by y_add
fig_pos(4) = fig_pos(4) - old_y_add + y_add; % make longer by y_add
set(edit_data.handles.edit_fig,'Position',fig_pos)

for j = 1:length(edit_data.handles.to_move)
    pos = get(edit_data.handles.to_move(j),'Position');
    pos(2) = pos(2) - old_y_add + y_add; % move down by y_add
    set(edit_data.handles.to_move(j),'Position',pos)
end

for i = 1:num_pts
    
    y = base_ht + (num_pts-i)*(edit_ht + space_ht);
    
    edit_data.handles.pt_text(i) = uicontrol(...
        'Style','text',...
        'Tag','pt_text',...
        'Units','characters',...
        'Position',[x_text y 12 1.5],...
        'BackgroundColor',b_color,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',10,...
        'FontWeight','normal',...
        'HorizontalAlignment','right',...
        'String',['Point ',num2str(i)]...
        );
    
    edit_data.handles.pt_edit(i) = uicontrol(...
        'Style','edit',...
        'Tag','pt_edit',...
        'Units','characters',...
        'Position',[x_edit y 20 1.75],...
        'BackgroundColor','white',...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',10,...
        'FontWeight','normal',...
        'HorizontalAlignment','left'...
        );
end

