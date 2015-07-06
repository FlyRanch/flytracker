% Function TYPE_MENU_CALLBACK
% 
% CALLING FUNCTION: type_menu_callback (Kine_v2_0)
% ACTIONS: Type menu (in edit dialog) callback for type 'MODEL'
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 17, 2004 by gwyneth

function type_menu_callback

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

%--------------------------------------------------------------------------
pict_height = 15;
pict_bottom_space = 3;
pict_top_space = 2;
pict_x = 7;
%--------------------------------------------------------------------------

edit_fig = findobj('Tag','edit_fig');
edit_data = guidata(edit_fig);

% Clear all objects specific to type, resize, set num_edit to 0   
to_clear = findobj('Parent',edit_fig);

for i = 1:length(edit_data.handles.to_keep)
    ind = find(to_clear ~= edit_data.handles.to_keep(i));
    to_clear = to_clear(ind);
end

delete(to_clear)

pos = get(edit_data.handles.num_edit,'Position');
y_sub = pos(2) - 6;

fig_pos = get(edit_data.handles.edit_fig,'Position');
fig_pos(2) = fig_pos(2) + y_sub; % move down by y_add
fig_pos(4) = fig_pos(4) - y_sub; % make longer by y_add
set(edit_data.handles.edit_fig,'Position',fig_pos)

for j = 1:length(edit_data.handles.to_move)
    pos = get(edit_data.handles.to_move(j),'Position');
    pos(2) = pos(2) - y_sub; % move down by y_add
    set(edit_data.handles.to_move(j),'Position',pos)
end

set(edit_data.handles.num_edit,'String','0')

% Get user to select model
cd([data.setup.mfile_path,filesep,'models'])
[model.name model.dir] = uigetfile('.mat','Choose a model to load');
model_name = model.name;
model_dir = model.dir;
if model_name == 0
    cd(data.setup.mfile_path)
    return
else
    cd(model_dir)
    load(model_name)
end
cd(data.setup.mfile_path)

if exist('object_type') == 0 % because not all models have this field yet
    warndlg('That is not a segment model, please choose again')
    return
end

if strcmp(object_type, 'segment_model') == 0
    warndlg('That is not a segment model, please choose again')
    return
end

num_pts = size(anchor_array,1);

% Make space for model picture
pos = get(edit_data.handles.num_edit,'Position');
y_add = pict_height;

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
y = y(2) + y(4) + pict_bottom_space;

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
        
%         % Plot parameter axes - NO PARAMETERS YET
%         for i = 1:size(param_array,1)
%             pts = param_array{i,2};
%             plot([pts(1) pts(3)],[pts(2) pts(4)],'g')
%             
%             param_x = (pts(1) + pts(3)) / 2;
%             param_y = (pts(2) + pts(4)) / 2;
%             plot_text(model_pict,param_x,param_y,0,param_array{i,1})
%         end
        
    case 3 %3D
        warndlg('Can''t handle 3D models yet')
        %plot3(coords(1,:),coords(2,:),coords(3,:))
        
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
          'String',['Selected Model: ',model_name])

guidata(controller, data)
edit_data.model = model;
guidata(edit_fig, edit_data)