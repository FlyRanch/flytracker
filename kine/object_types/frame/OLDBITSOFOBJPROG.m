%-------------------------------------------------------------------------


% % % Clear old drawn models --- edited by Qing Liu
% to_delete = findobj('Tag',cur_obj);
% delete(to_delete)
% 
% % Don't draw if all points not there
% undone = find(data.kine.(cur_obj).data.coords(:,:,frame) == 0);
% if isempty(undone)
% 
%     model_3d = [data.kine.(cur_obj).config.model_coords];
% 
% % Get spherical coordinates for selected model anchor points
% anch_1 = data.kine.(cur_obj).config.anchor_array{1,2};
% anch_2 = data.kine.(cur_obj).config.anchor_array{2,2};
% 
% mx1 = data.kine.(cur_obj).config.model_coords(1,anch_1);
% my1 = data.kine.(cur_obj).config.model_coords(2,anch_1);
% mz1 = data.kine.(cur_obj).config.model_coords(3,anch_1);
% mx2 = data.kine.(cur_obj).config.model_coords(1,anch_2);
% my2 = data.kine.(cur_obj).config.model_coords(2,anch_2);
% mz2 = data.kine.(cur_obj).config.model_coords(3,anch_2);
% 
% mdx = mx2 - mx1;   
% mdy = my2 - my1;
% mdz = mz2 - mz1;
% 
% [m_theta,m_phi,m_length] = cart2sph(mdx,mdy,mdz);
% 
% % Get spherical coordinates for digitized anchor points
% % cx1 = data.kine.(cur_obj).data.coords(1,1,frame);
% % cy1 = data.kine.(cur_obj).data.coords(1,2,frame);
% % cz1 = data.kine.(cur_obj).data.coords(1,3,frame);
% % cx2 = data.kine.(cur_obj).data.coords(2,1,frame);
% % cy2 = data.kine.(cur_obj).data.coords(2,2,frame);
% % cz2 = data.kine.(cur_obj).data.coords(2,3,frame);
% 
% % DEBUG: need to reverse definitions of c1 and c2 from previous models,
% % don't know why, but now it plots the frame in the correct orientation
% cx2 = data.kine.(cur_obj).data.coords(1,1,frame);
% cy2 = data.kine.(cur_obj).data.coords(1,2,frame);
% cz2 = data.kine.(cur_obj).data.coords(1,3,frame);
% cx1 = data.kine.(cur_obj).data.coords(2,1,frame);
% cy1 = data.kine.(cur_obj).data.coords(2,2,frame);
% cz1 = data.kine.(cur_obj).data.coords(2,3,frame);
% 
% cdx = cx2 - cx1;   
% cdy = cy2 - cy1;
% cdz = cz2 - cz1;
% 
% [c_theta,c_phi,c_length] = cart2sph(cdx,cdy,cdz);
% 
% % Scale factor is digitized units/model units
% s = c_length / m_length;    
% 
% %Get stored parameter value 
% alpha = (pi/180) * data.kine.(cur_obj).data.params(1,frame);      
% alpha = (2*pi)-alpha; %DEBUG NEW SYSTEM
% 
% % use rotation matrices in sequence to rotate points in wing model matrix
% %m_fit = rot3D(model_3d, -c_theta, -c_phi, -alpha); 
% %alpha = 270-alpha;
% m_fit = rot3D(model_3d, c_theta, -c_phi, alpha); %th=az=psi; phi=el=th; al=phi
% 
% % Scale wing to measured length (may be different lengths each frame)
% m_fit = s * m_fit;
% 
% % Translate wing model to right place
% m_fit(1,:) = m_fit(1,:) + cx1;
% m_fit(2,:) = m_fit(2,:) + cy1;
% m_fit(3,:) = m_fit(3,:) + cz1;
% 
% % use calibration data to calculate what the wing model should look like in each 2D window
% %%%DEBUG when add other calibration methods
% x = m_fit(1,:);
% y = m_fit(2,:);
% z = m_fit(3,:);
% 
% col_name = data.kine.(cur_obj).config.color;
% col_rgb = data.colors.(col_name);
% for i = 1:data.setup.cam_num
%     
%     [u(i).pts v(i).pts] = dlt_3D_to_2D( data.cal.coeff.(['DLT_',num2str(i)]), x, y, z );
%     %axes(dig_data.handles.(['cam',num2str(i)]));
%     cur_ax = dig_data.handles.(['cam',num2str(i)]);
%     %cam_pt(i) = plot(u(i).pts, v(i).pts, 'Tag', cur_obj);
%     cam_pt(i) = plot(cur_ax,u(i).pts(10:11), v(i).pts(10:11),'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible,'Tag', cur_obj);
%     other_pt(i) = plot(cur_ax,u(i).pts(14:15), v(i).pts(14:15),'Color',[.6 .6 .6],'Tag', cur_obj,'LineStyle',':');
%     other_pt(i) = plot(cur_ax,u(i).pts(10),v(i).pts(10),'Tag', cur_obj,'Marker','.','Color','m'); % LEFT wing hinge
%     
% end
% 
% % col_name = data.kine.(cur_obj).config.color;
% % col_rgb = data.colors.(col_name); % look color RGB values up in data.colors
% % set(cam_pt,'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible);
% 
% guidata(dig_fig,dig_data)
% 
% % dig_data = guidata(dig_fig);
% % dig_data.objects.(cur_obj).obj_program
% 
% % MODIFIED 01/16/07 -- to test quaternion rotation
% show_quat = 11;
% if show_quat == 1
%     p = [0 1 0 0]';       % starting position as a quaternion
%     [pf] = rotquat(p, c_theta, -c_phi, alpha);
%     plotquat(pf, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 2 % THIS ONE IS INCORRECT
%     p = [0 1 0 0]';       % starting position as a quaternion
%     [pf] = rotquat2(p, c_theta, -c_phi, alpha);
%     plotquat(pf, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 3
%     pf = rotquat_matrix(model_3d, c_theta, -c_phi, alpha);
%     plotquat_matrix(pf, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 4 % Use rotated euler angle scheme
%     addpath C:\MatlabRoot\Project_1\DigAnalysis_20070120\EulerAngleAnalysis20070907
%     p = [0 1 0 0]'; %not really needed
%     [pf, q] = rotquat(p, c_theta, -c_phi, alpha);
%     [mod_rot] = roteulMstar(model_3d,q);
%     plotquat_matrix(mod_rot, s, cur_obj, cx1, cy1, cz1)
%     
%     eul = quat2eul_vout(q);
%     [eul(3);-eul(2);360-eul(1)]
% elseif show_quat == 5 % Use rotated euler angle scheme
%     addpath C:\MatlabRoot\Project_1\DigAnalysis_20070120\EulerAngleAnalysis20070907
%     p = [0 1 0 0]'; %not really needed
%     [pf, q] = rotquat(p, c_theta, -c_phi, alpha);
%     [mod_rot] = roteulMstar5(model_3d,q);
%     plotquat_matrix(mod_rot, s, cur_obj, cx1, cy1, cz1)
%     
%     eul = quat2eul_vout(q);
%     [eul(3);-eul(2);360-eul(1)]
% elseif show_quat == 6 % Rotate using the input euler angles back calculated from quaternion
%     addpath C:\MatlabRoot\Project_1\DigAnalysis_20070120\EulerAngleAnalysis20070907
%     p = [0 1 0 0]'; %not really needed
%     [pf, q] = rotquat(p, c_theta, -c_phi, alpha);
%     [mod_rot] = roteulM(model_3d,q);
%     %[mod_rot] = rot3D(model_3d,eul(3)*pi/180,eul(2)*pi/180,eul(1)*pi/180);
%     plotquat_matrix(mod_rot, s, cur_obj, cx1, cy1, cz1)
%         eul = quat2eul_vout(q);
% 
%         phi = eul(1)*pi/180;
%         theta = eul(2)*pi/180;
%         psi = eul(3)*pi/180;   
%         [pf, q_test] = rotquat(p, psi, theta, phi+2*pi);
%     
%     [alpha*180/pi eul(1)+360;-c_phi*180/pi eul(2);c_theta*180/pi eul(3)]
%     [q q_test]
% elseif show_quat == 7 % Use rotated euler angle scheme
%     addpath C:\MatlabRoot\Project_1\DigAnalysis_20070120\EulerAngleAnalysis20070907
%     p = [0 1 0 0]'; %not really needed
%     [pf, q] = rotquat(p, c_theta, -c_phi, alpha);
%     [mod_rot] = roteulMstar3(model_3d,q);
%     plotquat_matrix(mod_rot, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 8 % Will's scheme
%     % Start with rotated frame and see if it gets to the correct place
%     addpath C:\MatlabRoot\Project_1\DigAnalysis_20070120\EulerAngleAnalysis20070907
%     p = [0 1 0 0]'; %not really needed
%     [pf, q] = rotquat(p, c_theta, -c_phi, alpha);
%     new_model = load('C:\MATLAB7\work\kine_v2_1\models\fly_frame_20061120_rotated');
%     plotquat_matrix(new_model.coords, s, cur_obj, cx1, cy1, cz1)
%     [mod_rot] = roteulMstar8(new_model.coords,q);
%     plotquat_matrix(mod_rot, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 9 % scheme with body axes defined properly for negative z
%     new_model = load('C:\MATLAB7\work\kine_v2_1\models\fly_frame_20061120_rotated_roll');
%     mod_rot = rot3D(new_model.coords, c_theta, -c_phi, pi-(2*pi-alpha));
%     alpha*180/pi
%     plotquat_matrix(mod_rot, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 10 % Quaternion with new alpha scheme
%     p = [0 1 0 0]';       % starting position as a quaternion
%     [pf] = rotquat(p, c_theta, -c_phi, pi-(2*pi-alpha));
%     plotquat(pf, s, cur_obj, cx1, cy1, cz1)
% elseif show_quat == 11 % Quaternion matrix with new alpha scheme
%     new_model = load('C:\MATLAB7\work\kine_v2_1\models\fly_frame_20061120_rotated_roll');
%     [pf, q] = rotquat_matrix(new_model.coords, c_theta, -c_phi, pi-(2*pi-alpha));
%     plotquat_matrix(pf, s, cur_obj, cx1, cy1, cz1)
%     data.test.q(:,frame) = q;
%     data.test.eul_xyz(:,frame) = [c_theta; -c_phi; pi-(2*pi-alpha)];
%     guidata(controller,data)
% end
% 
% 
% end
% 
% % Store match point values for the frame
% 
% figure(dig_fig)