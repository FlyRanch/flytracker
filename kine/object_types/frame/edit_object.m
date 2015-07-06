% Function EDIT_OBJECT
% 
% CALLING FUNCTION: Callback for Edit button in obj_fig for 'FRAME'
%                   objects
% ACTIONS: Opens edit object dialog then sets values of uicontrols
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 17, 2004 by gwyneth

function edit_object

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

addpath(data.setup.mfile_path)
%--------------------------------------------------------------------------
pict_height = 15;
pict_bottom_space = 3;
pict_top_space = 2;
pict_x = 7;

match_h = 2;
match_text_w = 10;
match_menu_w = 20;
match_above_space = 1;
match_below_space = 1;
match_btwn_space = 1;
match_horiz_space = 2;
az = -10.50;
el = 26;
%--------------------------------------------------------------------------

% Set value of uicontrols based on currently selected object
obj_fig = findobj('Tag','obj_fig');
edit_fig = findobj('Tag','edit_fig');

obj_data = guidata(obj_fig);
edit_data = guidata(edit_fig);

obj_name = get(edit_data.handles.name_edit,'String');

% Get user to select model
model_name = obj_data.kine.(obj_name).config.model_name;
coords = obj_data.kine.(obj_name).config.model_coords;
anchor_array = obj_data.kine.(obj_name).config.anchor_array;
param_array = obj_data.kine.(obj_name).config.param_array;
match_points = obj_data.kine.(obj_name).config.match_points;
match_array = obj_data.kine.(obj_name).config.match_array;

num_pts = size(anchor_array,1);

% Make space for model picture & match point assignments
pos = get(edit_data.handles.num_edit,'Position');

num_match_pts = length(match_points);
m_add = num_match_pts*(match_h + match_btwn_space) + match_above_space + match_below_space;
y_add = pict_height + m_add;

fig_pos = get(edit_data.handles.edit_fig,'Position');
fig_pos(2) = fig_pos(2) - y_add; % move down by y_add
fig_pos(4) = fig_pos(4) + y_add; % make longer by y_add
set(edit_data.handles.edit_fig,'Position',fig_pos)

for j = 1:length(edit_data.handles.to_move)
    pos = get(edit_data.handles.to_move(j),'Position');
    pos(2) = pos(2) + y_add; % move down by y_add
    set(edit_data.handles.to_move(j),'Position',pos)
end

y = get(edit_data.handles.ok_edit_button,'Position');
y = y(2) + y(4) + pict_bottom_space + m_add;

h = pict_height - pict_bottom_space - pict_top_space;
w = fig_pos(3) - 2*pict_x;

% Create axes and plot model with anchor/param labels
model_pict = axes('Units','characters','Position',[pict_x y w h]);
set(model_pict,'FontSize',8)
hold on

switch size(coords,1)
    case 2 %2D
        % Plot model coords
        plot(coords(1,:),coords(2,:));
        
        % Plot anchor points and label
        for i = 1:size(anchor_array,1)
            plot(coords(1,anchor_array{i,2}),coords(2,anchor_array{i,2}),'ro')
            plot_text(model_pict,coords(1,anchor_array{i,2}),coords(2,anchor_array{i,2}),5,anchor_array{i,1})
        end
        
        % Plot parameter axes
        for i = 1:size(param_array,1)
            pts = param_array{i,2};
            plot([pts(1) pts(3)],[pts(2) pts(4)],'g')
            
            param_x = (pts(1) + pts(3)) / 2;
            param_y = (pts(2) + pts(4)) / 2;
            plot_text(model_pict,param_x,param_y,0,param_array{i,1})
        end
        
    case 3 %3D
        plot3(coords(1,:),coords(2,:),coords(3,:))
        axis equal
        view(az,el)
        
        % Plot anchor points and label
        for i = 1:size(anchor_array,1)
            plot3(coords(1,anchor_array{i,2}),coords(2,anchor_array{i,2}),coords(3,anchor_array{i,2}),'ro')
            plot_text3(model_pict,coords(1,anchor_array{i,2}),coords(2,anchor_array{i,2}),coords(3,anchor_array{i,2}),anchor_array{i,1},25,'-z')
        end

        % Plot match points and label
        for i = 1:length(match_points)
            plot3(coords(1,match_points(i)),coords(2,match_points(i)),coords(3,match_points(i)),'bo')
            plot_text3(model_pict,coords(1,match_points(i)),coords(2,match_points(i)),coords(3,match_points(i)),char(i+96),5,'-x')
        end

        % Plot parameter axes
        for i = 1:size(param_array,1)
            pts = param_array{i,2};
            plot3(pts(:,1),pts(:,2),pts(:,3),'g')
            
            param_x = mean(pts(:,1));
            param_y = mean(pts(:,2));
            param_z = mean(pts(:,3));
            plot_text3(model_pict,param_x,param_y,param_z,param_array{i,1},25,'+z')
        end

        
    otherwise
        warndlg('uh-oh!')
        return
end

% Set num_edit
set(edit_data.handles.num_edit,'String',num2str(num_pts),'Style','text','HorizontalAlignment','center')

% Display model name
uicontrol('Style','text',...
          'Units','characters',...
          'Position',[0 y+h+1 fig_pos(3) 2],...
          'BackgroundColor',get(edit_data.handles.edit_fig,'Color'),...
          'FontSize',10,...
          'Tag','model name text',...
          'String',['Selected Frame: ',model_name])
      
% Display match choices
y = get(edit_data.handles.ok_edit_button,'Position');
y = y(2) + y(4) + m_add - match_above_space;

objects = fieldnames(obj_data.kine);

for i = 1:num_match_pts
    
    my = y - i*(match_h+match_btwn_space) + match_btwn_space;
    cur_obj = match_array{i,2};
    obj_val = strmatch(match_array{i,2},objects);
    points = obj_data.kine.(cur_obj).config.points;
    pt_val = strmatch(match_array{i,3},points);
    
    uicontrol('Style','text',...
              'Units','characters',...
              'Position',[pict_x my match_text_w match_h],...
              'BackgroundColor',get(edit_data.handles.edit_fig,'Color'),...
              'FontSize',10,...
              'Tag','match name text',...
              'String',[char(i+96),' )'])
          
    uicontrol('Style','popupmenu',...
              'Units','characters',...
              'Position',[(pict_x+match_text_w+match_horiz_space) my match_menu_w match_h],...
              'BackgroundColor','white',...
              'FontSize',10,...
              'Tag',['ob_menu_',num2str(i)],...
              'String',objects,...
              'Value',obj_val,...
              'UserData',{'match_menu_callback'},...
              'Callback','edit_fig_type_specific_callback')
          
    uicontrol('Style','popupmenu',...
              'Units','characters',...
              'Position',[(pict_x+match_text_w+match_menu_w+2*match_horiz_space) my match_menu_w match_h],...
              'BackgroundColor','white',...
              'FontSize',10,...
              'Tag',['pt_menu_',num2str(i)],...
              'String',points,...
              'Value',pt_val)
          
end

guidata(controller,data);

edit_data.model_name = model_name;

guidata(edit_fig, edit_data);