% Function OBJ_PROGRAM(cur_obj)
% 
% CALLING FUNCTION: model_params, param_edit_callback, param_button_callback
% ACTIONS: Draw model on the camera views with current parameter settings
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: November 16, 2004 by gwyneth

function obj_program(cur_obj)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

cur_type = data.kine.(cur_obj).config.type;
obj_num = get(data.handles.cur_object_menu,'Value');
frame = str2num(get(data.handles.frame_box,'String'));
num = size(data.kine.(cur_obj).config.param_array,1);%%%DEBUG for now this is hard-wired for 1 parameter

% % Clear old drawn models --- edited by Qing Liu
to_delete = findobj('Tag',cur_obj);
delete(to_delete)

% Don't draw if all points not there
undone = find(data.kine.(cur_obj).data.coords(:,:,frame) == 0);
if isempty(undone)

    model_3d = [data.kine.(cur_obj).config.model_coords; zeros(1,length(data.kine.(cur_obj).config.model_coords)) ];

% Bend the model frame, scale (For double bend, bend axes are hard-wired
% into bending program)
theta = strmatch('theta',data.kine.(cur_obj).config.param_array(:,1),'exact');
theta2 = strmatch('theta2',data.kine.(cur_obj).config.param_array(:,1),'exact');
% dis = strmatch('d',data.kine.(cur_obj).config.param_array(:,1),'exact');
% ang = strmatch('a',data.kine.(cur_obj).config.param_array(:,1),'exact');

% d = data.kine.(cur_obj).data.params(dis,frame);%0.15; % Could make this a parameter
% a = data.kine.(cur_obj).data.params(ang,frame);%70;  % Could make this a parameter
th1 = data.kine.(cur_obj).data.params(theta,frame);
th2 = data.kine.(cur_obj).data.params(theta2,frame);

[wing_bent, ax1, ax2] = wing_double_bend(model_3d, th1, th2);

% Get spherical coordinates for selected model anchor points
anch_1 = data.kine.(cur_obj).config.anchor_array{1,2};
anch_2 = data.kine.(cur_obj).config.anchor_array{2,2};

mx1 = data.kine.(cur_obj).config.model_coords(1,anch_1);
my1 = data.kine.(cur_obj).config.model_coords(2,anch_1);
mx2 = data.kine.(cur_obj).config.model_coords(1,anch_2);
my2 = data.kine.(cur_obj).config.model_coords(2,anch_2);

mdx = mx2 - mx1;   
mdy = my2 - my1;
mdz = 0; %models must be planar

[m_theta,m_phi,m_length] = cart2sph(mdx,mdy,mdz); %(only use m_length)

% Get spherical coordinates for digitized anchor points
cx1 = data.kine.(cur_obj).data.coords(1,1,frame);
cy1 = data.kine.(cur_obj).data.coords(1,2,frame);
cz1 = data.kine.(cur_obj).data.coords(1,3,frame);
cx2 = data.kine.(cur_obj).data.coords(2,1,frame);
cy2 = data.kine.(cur_obj).data.coords(2,2,frame);
cz2 = data.kine.(cur_obj).data.coords(2,3,frame);

cdx = cx2 - cx1;   
cdy = cy2 - cy1;
cdz = cz2 - cz1;

[c_theta,c_phi,c_length] = cart2sph(cdx,cdy,cdz);

% Scale factor is digitized units/model units
s = c_length / m_length;    

%Get stored parameter value 
alpha = (pi/180) * data.kine.(cur_obj).data.params(1,frame);      

% use rotation matrices in sequence to rotate points in wing model matrix
m_fit = rot3D(wing_bent, -c_theta, -c_phi, -alpha);    

% Scale wing to measured length (may be different lengths each frame)
m_fit = s * m_fit;

% Translate wing model to right place
m_fit(1,:) = m_fit(1,:) + cx1;
m_fit(2,:) = m_fit(2,:) + cy1;
m_fit(3,:) = m_fit(3,:) + cz1;

% use calibration data to calculate what the wing model should look like in each 2D window
%%%DEBUG when add other calibration methods
x = m_fit(1,:);
y = m_fit(2,:);
z = m_fit(3,:);

[ u_1, v_1 ] = dlt_3D_to_2D( data.cal.coeff.DLT_1, x, y, z );
[ u_2, v_2 ] = dlt_3D_to_2D( data.cal.coeff.DLT_2, x, y, z );
[ u_3, v_3 ] = dlt_3D_to_2D( data.cal.coeff.DLT_3, x, y, z );

axes(dig_data.handles.cam1)
cam1_pt = plot(u_1,v_1,'Tag',cur_obj);

axes(dig_data.handles.cam2)
cam2_pt = plot(u_2,v_2,'Tag',cur_obj);

axes(dig_data.handles.cam3)
cam3_pt = plot(u_3,v_3,'Tag',cur_obj);

% Draw bending axes
axes(dig_data.handles.cam1)
cam1_ax1 = plot(u_1(ax1),v_1(ax1),':','Tag',cur_obj);
cam1_ax2 = plot(u_1(ax2),v_1(ax2),':','Tag',cur_obj);

axes(dig_data.handles.cam2)
cam2_ax1 = plot(u_2(ax1),v_2(ax1),':','Tag',cur_obj);
cam2_ax2 = plot(u_2(ax2),v_2(ax2),':','Tag',cur_obj);

axes(dig_data.handles.cam3)
cam3_ax1 = plot(u_3(ax1),v_3(ax1),':','Tag',cur_obj);
cam3_ax2 = plot(u_3(ax2),v_3(ax2),':','Tag',cur_obj);

% Save handles of ploted marks into dig_fig.objects structure
% dig_data.objects.(cur_obj).obj_program(1,1:3) = [cam1_pt cam2_pt cam3_pt];
    
col_name = data.kine.(cur_obj).config.color;
col_rgb = data.colors.(col_name); % look color RGB values up in data.colors
set([cam1_pt cam2_pt cam3_pt],'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible)
set([cam1_ax1 cam2_ax1 cam3_ax1 cam1_ax2 cam2_ax2 cam3_ax2],'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible)

guidata(dig_fig,dig_data)

% dig_data = guidata(dig_fig);
% dig_data.objects.(cur_obj).obj_program
end
