% Function EDIT_OBJECT
% 
% CALLING FUNCTION: Callback for Edit button in obj_fig for 'MODEL'
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
% param_array = obj_data.kine.(obj_name).config.param_array; NO PARAMETERS YET

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

guidata(controller,data);

edit_data.model_name = model_name;

guidata(edit_fig, edit_data);