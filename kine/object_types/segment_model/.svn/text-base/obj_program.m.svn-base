% Function OBJ_PROGRAM(cur_obj)
% 
% CALLING FUNCTION: model_params, param_edit_callback, param_button_callback
% ACTIONS: Draw model on the camera views with current parameter settings
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: July 28, 2004 by gwyneth

function obj_program(cur_obj)

%--------------------------------------------------------------------------
%                      FEMUR         TIBIA      TARSAL SEGMENTS
% (1)trochanter  o--------------o------------o------------------o (4)tarsus
%                   femur-tibia(2)          (3)tibia-tarsus

tro = 1;
ft = 2;
tt = 3;
tar = 4;
x = 1;
y = 2;
z = 3;
%--------------------------------------------------------------------------

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

cur_type = data.kine.(cur_obj).config.type;
obj_color = data.kine.(cur_obj).config.color;
obj_num = get(data.handles.cur_object_menu,'Value');
frame = str2num(get(data.handles.frame_box,'String'));
% num = size(data.kine.(cur_obj).config.param_array,1);NO PARAMETERS YET
coords = data.kine.(cur_obj).data.coords(:,:,frame);
model_coords = data.kine.(cur_obj).config.model_coords;
num_cams = data.setup.cam_num;

% Clear old drawn models --- edited by Qing Liu
to_delete = findobj('Tag',cur_obj);
delete(to_delete)

% DRAW TIBIA
istibia = 'no';
if coords(ft,x) ~= 0 & coords(tt,x) ~= 0
    l_tibia = pt_dist3( coords(ft,x), coords(ft,y), coords(ft,z), coords(tt,x), coords(tt,y), coords(tt,z) );
    drawseg(cur_obj,frame,ft:tt,num_cams)
    istibia = 'yes';
end    
    
% DRAW FEMUR
if strcmp(istibia,'yes') & coords(tro,x) ~= 0
    t = ( model_coords(1,2) / (model_coords(1,3)-model_coords(1,2)) ) * l_tibia;     % scaled length of femur
    data.kine.(cur_obj).data.coords(tro,:,frame) = projpt( coords(ft,:), coords(tro,:), t );
    guidata(controller,data)
    drawseg(cur_obj,frame,tro:ft,num_cams)
elseif ~strcmp(istibia,'yes') & coords(tro,x) ~= 0
    set(dig_data.objects.(cur_obj).points(tro,:),'MarkerEdgeColor','w')
end

% DRAW TARSAL SEGMENT(S)
if strcmp(istibia,'yes') & coords(tar,x) ~= 0
    t = ( (1-model_coords(1,3)) / (model_coords(1,3)-model_coords(1,2)) ) * l_tibia; % scaled length of tarsals
    data.kine.(cur_obj).data.coords(tar,:,frame) = projpt( coords(tt,:), coords(tar,:), t );
    guidata(controller,data)
    drawseg(cur_obj,frame,tt:tar,num_cams)
elseif ~strcmp(istibia,'yes') & coords(tar,x) ~= 0
    set(dig_data.objects.(cur_obj).points(tar,:),'MarkerEdgeColor','w')
end

%--------------------------------------------------------------------------
function drawseg(cur_obj,frame,pts,num_cams)

controller = findobj('Tag','controller'); data = guidata(controller);
dig_fig = findobj('Tag','dig_fig'); dig_data = guidata(dig_fig);

col_name = data.kine.(cur_obj).config.color;
col_rgb = data.colors.(col_name); % look color RGB values up in data.colors

for cam = 1:num_cams
    [ u, v ] = dlt_3D_to_2D( data.cal.coeff.(['DLT_',num2str(cam)]),...
        data.kine.(cur_obj).data.coords(pts,1,frame),...
        data.kine.(cur_obj).data.coords(pts,2,frame),...
        data.kine.(cur_obj).data.coords(pts,3,frame) );
    
    %axes( dig_data.handles.(['cam',num2str(cam)]) )
    cur_ax = dig_data.handles.(['cam',num2str(cam)]);
    plot( cur_ax, u, v, 'Tag', cur_obj, 'Color', col_rgb, 'Visible',data.kine.(cur_obj).config.visible )
end

%--------------------------------------------------------------------------
function endpt = projpt( anchor, segpt, t )

x = 1; y = 2; z = 3;

l_seg = pt_dist3( anchor(x), anchor(y), anchor(z), segpt(x), segpt(y), segpt(z) );
v_seg = [ segpt(x)-anchor(x) segpt(y)-anchor(y) segpt(z)-anchor(z) ] / l_seg;

endpt = anchor + t*v_seg;
%--------------------------------------------------------------------------