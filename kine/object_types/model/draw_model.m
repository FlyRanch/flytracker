% Function DRAW_MODEL
% 
% CALLING FUNCTION: obj_function
% ACTIONS: Draw model on the camera views with current parameter settings
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: October 3, 2007 by gwyneth

function draw_model(cur_obj)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

frame = str2double(get(data.handles.frame_box,'String'));

% Clear old drawn models
to_delete = findobj('Tag',cur_obj);
delete(to_delete)

% Get rotated & translated model points to plot from Kine
m_fit = data.kine.(cur_obj).data.model_coords(:,:,frame);

% Use calibration data project model into each 2D window (DEBUG when add
% other calibration methods)
x = m_fit(1,:);
y = m_fit(2,:);
z = m_fit(3,:);

col_name = data.kine.(cur_obj).config.color;
col_rgb = data.colors.(col_name);
for i = 1:data.setup.cam_num
    
    [u v] = dlt_3D_to_2D( data.cal.coeff.(['DLT_',num2str(i)]), x, y, z );
    cur_ax = dig_data.handles.(['cam',num2str(i)]);
    plot(cur_ax,u, v,'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible,'Tag', cur_obj);
    %plot(cur_ax,u(14:15), v(14:15),'Color',[.6 .6 .6],'Tag', cur_obj,'LineStyle',':');
    %plot(cur_ax,u(10),v(10),'Tag', cur_obj,'Marker','.','Color','m'); % LEFT wing hinge
    
end

guidata(dig_fig,dig_data)
figure(dig_fig)