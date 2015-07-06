% Function DRAW_BCVIEW(cur_obj)
% 
% CALLING FUNCTION: draw_it
% ACTIONS: Plots an object in body-centered view
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: September 27, 2007 by gwyneth

function draw_bcview(cur_obj, frame)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Get object to body center around
bc_obj = data.bc.bc_obj;

% Get the points to rotate
draw_model = 0;
if isfield(data.kine.(cur_obj).data,'model_coords')
    m_fit = data.kine.(cur_obj).data.model_coords(:,:,frame);
    draw_model = 1;
end
pts = data.kine.(cur_obj).data.coords(:,:,frame)';

% Get the translation vector
v_trans = data.kine.(bc_obj).data.v_trans(:,frame);

% Get the quaternion to use for rotation
q = data.kine.(bc_obj).data.quat(:,frame);

% Translate the points
pts = pts - repmat(v_trans,1,size(pts,2));
if draw_model == 1;
    m_fit = m_fit - repmat(v_trans,1,size(m_fit,2));
end

% Rotate the points!
p_bc = quatRot_lab2bod(q,pts);
if draw_model == 1;
    m_bc = quatRot_lab2bod(q,m_fit);
end

% Plot in 3D, but then restrict axes to specified views
side = dig_data.handles.body_plot_side;
rear = dig_data.handles.body_plot_rear;
top  = dig_data.handles.body_plot_top;
h = [side, rear, top];

xmin = -2;
xmax = 2;
ymin = xmin;
ymax = xmax;
zmin = xmin;
zmax = xmax;
h_views = [180 0; -90 0; 180 90];
xtxt = {'x','y','x'};
ytxt = {'z','z','y'};

col_name = data.kine.(cur_obj).config.color;
col_rgb = data.colors.(col_name);

for i = 1:length(h)

    for j = 1:size(p_bc,2)
        if p_bc(1,j)~=0 & p_bc(2,j)~=0
            plot3(h(i), p_bc(1,:), p_bc(2,:), p_bc(3,:),...
                'Marker', 'o', 'MarkerSize', 2, 'LineStyle','none','Color', col_rgb,...
                'Visible',data.kine.(cur_obj).config.visible,'Tag', cur_obj);
        end
    end
    
    if draw_model == 1
        plot3(h(i), m_bc(1,:), m_bc(2,:), m_bc(3,:),'Color',col_rgb,...
            'Visible',data.kine.(cur_obj).config.visible,'Tag', cur_obj);
    end
    %set(h(i),'zdir','reverse','xdir','reverse','ydir','reverse')
    
    axis(h(i),[xmin xmax ymin ymax zmin zmax])
    view(h(i),h_views(i,1),h_views(i,2))
    
%     xlabel(h(i),xtxt{i})
%     ylabel(h(i),ytxt{i})
end
