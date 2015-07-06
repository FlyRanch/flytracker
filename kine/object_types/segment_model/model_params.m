% Function MODEL_PARAMS
% 
% CALLING FUNCTION: dig_setup (model), fit_model button callback
% ACTIONS: Displays model parameters and fits model
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 30, 2004 by gwyneth

function model_params

controller = findobj('Tag','controller');
data = guidata(controller);

% Check to see if all the anchor points of the current object are digitized
% for the current frame
cur_obj = get_string(data.handles.cur_object_menu);
frame = str2num(get(data.handles.frame_box,'String'));

undone = find(data.kine.(cur_obj).data.coords(:,:,frame) == 0);
if isempty(undone) == 0
    
    warndlg('The model cannot be fit until all the anchor points are dtermined for this frame')
    return
    
end

% Delete the fit model button
call_param_button = findobj('Tag','call_param_button');
delete(call_param_button)

%%%DEBUG
to_delete = findobj('Tag','obj_control');
delete(to_delete)

% Format the parameter controls according to how many parameters there are
num_params = size(data.kine.(cur_obj).config.param_array,1);

    % NOTE: can't just set up the obj controls in character units because
    % main windows will be a different size on every screen, so need to
    % find the dimensions of the the object controls frame (obj_frame) and 
    % use that to base the location of the object controls...and the size, 
    % really...
    
obj_frame_pos = get(findobj('Tag','obj_frame'),'Position');
xmin = obj_frame_pos(1);
xmax = obj_frame_pos(1) + obj_frame_pos(3);
ymin = obj_frame_pos(2);
ymax = obj_frame_pos(2) + obj_frame_pos(4);

horiz_marg = xmin + 0.03;
vert_marg = ymin + 0.0075;
horiz_btwn = 0.03;
vert_btwn = 0.0075;
l_edit = 0.09;
h_edit = 0.025;
l_slider = 2*l_edit;
h_slider = h_edit;
text_gap = 0.02;
l_text = 0.10;
h_text = h_edit;

for i = 1:num_params
    
    if i <= 5
        x = horiz_marg;
        y = vert_marg + (5 - i)*(h_edit + vert_btwn);
        
    elseif i > 5 & i <= 10
        x = horiz_marg + 1*(l_slider + l_edit + text_gap + l_text + horiz_btwn);
        y = vert_marg + (10 - i)*(h_edit + vert_btwn);
        
    elseif i > 10 
        
        warndlg('Can only show controls for the first 10 parameters')
        
%         & i <= 15
%         x = horiz_marg + 2*(l_slider + l_edit + text_gap + l_text + horiz_btwn);
%         y = vert_marg + (15 - i)*(h_edit + vert_btwn);
        
    end
    
    figure(controller)
    
    s_max = data.kine.(cur_obj).config.param_array{i,3}(2);
    s_min = data.kine.(cur_obj).config.param_array{i,3}(1);
    s_step(1) = 1/(s_max - s_min );
    s_step(2) = 5/(s_max - s_min );
    
    uicontrol(...
        'Style','slider',...
        'Tag','obj_control',...
        'Units','normalized',...
        'Position',[x y l_slider h_slider],...
        'BackgroundColor',data.colors.background,...
        'Max',s_max,...
        'Min',s_min,....
        'SliderStep',s_step,...
        'Value',data.kine.(cur_obj).data.params(i,frame),...
        'UserData',i,...
        'Callback','param_slider_callback'...
        );

    uicontrol(...
        'Style','edit',...
        'Tag','obj_control',...
        'Units','normalized',...
        'Position',[x+l_slider y l_edit h_edit],...
        'BackgroundColor','white',...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',12,...
        'FontWeight','normal',...
        'String',round(data.kine.(cur_obj).data.params(i,frame)),...
        'UserData',i,...
        'Callback','param_edit_callback'...
        );

    uicontrol(...
        'Style','text',...
        'Tag','obj_control',...
        'Units','normalized',...
        'Position',[x+l_slider+l_edit+text_gap y-0.005 l_text h_text],...
        'BackgroundColor',data.colors.background,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',11,...
        'FontWeight','normal',...
        'HorizontalAlignment','left',...
        'UserData',i,...
        'String',data.kine.(cur_obj).config.param_array{i,1}....
        );    
    
end

obj_program(cur_obj)  %%DEBUG: path?