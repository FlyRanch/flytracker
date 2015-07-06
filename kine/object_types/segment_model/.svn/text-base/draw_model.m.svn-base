% OBSOLETE, see obj_program

% Function DRAW_MODEL
% 
% CALLING FUNCTION: param_edit_callback, param_button_callback
% ACTIONS: Draw model on the camera views with current parameter settings
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: June 1, 2004 by gwyneth

function draw_model(cur_obj)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

%cur_obj = get_string(data.handles.cur_object_menu);
cur_type = data.kine.(cur_obj).config.type;
% obj_num = get(data.handles.cur_object_menu,'Value');
    objs = fieldnames(data.kine);
    obj_num = strmatch(cur_obj,objs);

frame = str2num(get(data.handles.frame_box,'String'));
num = size(data.kine.(cur_obj).config.param_array,1);%%%DEBUG for now this is hard-wired for 1 parameter

model_3d = [data.kine.(cur_obj).config.model_coords; zeros(1,length(data.kine.(cur_obj).config.model_coords)) ];

% Clear old drawn models --- edited by Qing Liu
if isfield(dig_data,'objects')
    if size(dig_data.objects,2) >= obj_num
        if isfield(dig_data.objects(obj_num),cur_type)
            delete(dig_data.objects(obj_num).(cur_type)(:))
        end
    end
end

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

[m_theta,m_phi,m_length] = cart2sph(mdx,mdy,mdz);

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
m_fit = rot3D(model_3d, -c_theta, -c_phi, -alpha);    

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
cam1_pt = plot(u_1,v_1);

axes(dig_data.handles.cam2)
cam2_pt = plot(u_2,v_2);

axes(dig_data.handles.cam3)
cam3_pt = plot(u_3,v_3);

% Save handles of ploted marks into dig_fig.objects structure
dig_data.objects(obj_num).(cur_type)(1,1:3) = [cam1_pt cam2_pt cam3_pt];
col_name = data.kine.(cur_obj).config.color;
col_rgb = data.colors.(col_name); % look color RGB values up in data.colors
set([cam1_pt cam2_pt cam3_pt],'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible)

guidata(dig_fig,dig_data)

disp('A wing model was drawn')