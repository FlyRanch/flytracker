% view saved kinematics that have been tracked manually


%clear 
%close all
fignum = 3;
plotscale = 1;

addpath kine/kine_math


PAR.dt = 1/6000;
PAR.etamax = 0;
PAR.c = 4;
PAR.pixpermm = 1;
PAR.numfly = 1;
PAR.L1 = 25; %# of steps for body along length
PAR.L2 = 10; %# of steps for head along length
PAR.L3 = 30; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

PAR.numcam = 3;
PAR.videopath = 'video/';
PAR.filetag = 'exp101';

%load kine/SavedKinematics/Wings/exp083_forEbraheem_corAlpha_splined_adj.mat
load kine/SavedKinematics/Wings/exp101_ebraheem.mat

% Get Calibration
for j = 1:PAR.numcam
    PAR.DLT(:,j) = data.cal.coeff.(['DLT_' num2str(j)]);
end

load flygenmod
params = flygenmod;

PAR.statedim = 12;

State = zeros(PAR.statedim,data.images.frames);
%frames = 400:490;
for k = 1:258%data.images.frames
    %===============================================================
%% Get the body orientation

    %Taken from Gwyneth's 'obj_program.m'
    
    % Use when data saved in old format
%     alpha = (pi/180) * data.kine.body.data.params(1,k);  
%     %alpha = pi - alpha;
%     
%     c1 = data.kine.body.data.coords(1,:,k);
%     c2 = data.kine.body.data.coords(2,:,k);
%     dc = c1 - c2;
%     [psi,theta,c_length] = cart2sph(dc(1),dc(2),dc(3));
%     
%     alpha = alpha;%pi-alpha;%(2*pi)-alpha; %DEBUG NEW SYSTEM: for new will just be alpha
%     theta = -theta; % reverse theta because cart2sph defines elevation above the xy-plane as positive, but this is a negative rotation around the y-axis
%     psi = psi; % don't reverse because rotate from body to lab frame
% 
%     q_body = eulzyx2quat(alpha,theta,psi);
    
    q_body = data.kine.body.data.quat(:,k);
    
    %I store the scalar part of quaternion at the end
    q_body = [q_body(2:4); q_body(1)];
    
    % add a rotation by pi along the roll axis since our model has the body
    % fixed frame with z-axis pointing ventral.
    q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);
    
%     % I want translation from body to fixed
%     q_body = [-q_body(1:3); q_body(4)];
    
    % Get Translation
    % Use this code when saved in old format
%     model_3d = data.kine.body.config.model_coords;
%     
%     % Get Model length (really this should always be 1)
%     if isfield(data.kine.body.config,'model_length')
%         m_length = data.kine.body.config.model_length;
%     else
%         anch_1 = data.kine.body.config.anchor_array{1,2};
%         anch_2 = data.kine.body.config.anchor_array{2,2};
% 
%         m1 = data.kine.body.config.model_coords(:,anch_1);
%         m2 = data.kine.body.config.model_coords(:,anch_2);
%         dm = m2 - m1;
% 
%         [m_theta,m_phi,m_length] = cart2sph(dm(1),dm(2),dm(3));
%         
%     end
%     % Scale factor is digitized units/model units
%     s = c_length / m_length;
% 
%     m_fit = quatRot_bod2lab(q_body, model_3d);
%     m_fit = scale3D(m_fit,s);
% 
%     pt_name = data.kine.body.config.anchor_array{1,1}; % use first anchor point
%     pt_col = data.kine.body.config.anchor_array{1,2};
%     kine_row = strmatch(pt_name, data.kine.body.config.points);
%     vb_anch = m_fit(:,pt_col);
%     vf_anch = data.kine.body.data.coords(kine_row,:,k);
% 
%     [m_fit, v_trans] = trans3D(m_fit, vb_anch, vf_anch');
%     
%     
%     T_body = v_trans;
    
    T_body = data.kine.body.data.v_trans(:,k);
    %Calculation translation vector -q_body*T10*q_body` + T_body
    %T10 is the approximate location of tailpt in our generative model
    %It is the next to last cross-section of model at the dorsal edge
    [xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(params.bodyctr,params.bodyrad,params.bodylen,PAR);

    %Say T10 is right on the tip of the thorax
    T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
    T10 = -T10;
    [T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),T_body,q_body);

    BodyTrans(:,k) = T_body;
    
%===============================================================
% Get the wing orientations
%% Left

% Use this section when the data is saved in the old format     
%     alpha = (pi/180) * data.kine.l_wing.data.params(1,k);%Get stored parameter value       
%     alpha = pi - alpha;
%     
%     tip = 2;
%     hinge = 1;
% 
%     c1 = data.kine.l_wing.data.coords(tip,:,k);
%     c2 = data.kine.l_wing.data.coords(hinge,:,k);
%     dc = c1 - c2;
%     [psi,theta,c_length] = cart2sph(dc(1),dc(2),dc(3));
%     
%     alpha = alpha;%pi-alpha;%(2*pi)-alpha; %DEBUG NEW SYSTEM: for new will just be alpha
%     theta = -theta; % reverse theta because cart2sph defines elevation above the xy-plane as positive, but this is a negative rotation around the y-axis
%     psi = psi; % don't reverse because rotate from body to lab frame
% 
%     q_Lwing = eulzyx2quat(alpha,theta,psi);
    
    %q_Lwing = data.kine.l_wing.data.quat(:,k);
    q_Lwing = data.kine.left_wing.data.quat(:,k);
    %I store the scalar part of quaternion at the end
    q_Lwing = [q_Lwing(2:4); q_Lwing(1)];

    % Premultiply rotation by this alignment quaternion that takes into account
    % the orientation of coordinate axis fixed to the left wing.
    q_Lwingaxisalign = quat2matNEW([0 -1 0;1 0 0;0 0 1]);
    q_Lwing = quatprod(q_Lwing,q_Lwingaxisalign);

%     % I want translation from body to fixed
%     q_Lwing = [-q_Lwing(1:3); q_Lwing(4)];

    % I will calculate the relative rotations
    % for the wings by multiplying the quaternions directly.
    % Just take the orientation because I assume that the wing is fixed at the
    % joint.
    q_Lwing_rel = quatprod([-q_body(1:3) ; q_body(4)],q_Lwing);
    
    T_L = data.kine.left_wing.data.v_trans(:,k);
    T_L = T_L(:);
   
    % Relative translation from body fixed axis to wing joint
    params.T_Lwing_rel_All(:,k) = qxform([-q_body(1:3) ; q_body(4)],T_L - T_body);

%% Right Wing

% Use this when data saved in old format
%     alpha = (pi/180) * data.kine.r_wing.data.params(1,k);%Get stored parameter value       
%     alpha = pi - alpha;
%     
%     tip = 2;
%     hinge = 1;
% 
%     c1 = data.kine.r_wing.data.coords(tip,:,k);
%     c2 = data.kine.r_wing.data.coords(hinge,:,k);
%     dc = c1 - c2;
%     [psi,theta,c_length] = cart2sph(dc(1),dc(2),dc(3));
%     
%     alpha = alpha;%pi-alpha;%(2*pi)-alpha; %DEBUG NEW SYSTEM: for new will just be alpha
%     theta = -theta; % reverse theta because cart2sph defines elevation above the xy-plane as positive, but this is a negative rotation around the y-axis
%     psi = psi; % don't reverse because rotate from body to lab frame
% 
%     q_Rwing = eulzyx2quat(alpha,theta,psi);

    %q_Rwing = data.kine.r_wing.data.quat(:,k);
    q_Rwing = data.kine.right_wing.data.quat(:,k);
    %I store the scalar part of quaternion at the end
    q_Rwing = [q_Rwing(2:4); q_Rwing(1)];

    % Premultiply rotation by this alignment quaternion that takes into account
    % the orientation of coordinate axis fixed to the right wing.
    q_Rwingaxisalign = quat2matNEW([0 1 0;1 0 0;0 0 -1]);
    q_Rwing = quatprod(q_Rwing,q_Rwingaxisalign);

%     % I want translation from body to fixed
%     q_Rwing = [-q_Rwing(1:3); q_Rwing(4)];

    
    %Just take the orientation because I assume that the wing is fixed at the
    %joint.
    q_Rwing_rel = quatprod([-q_body(1:3); q_body(4)],q_Rwing);

    T_R = data.kine.right_wing.data.v_trans(:,k);
    T_R = T_R(:);

    % Relative translation from body fixed axis to wing joint
    params.T_Rwing_rel_All(:,k) = qxform([-q_body(1:3) ; q_body(4)],T_R - T_body);

%% Calculate the twist representation for these quaternions
    R_body = quat2matNEW(q_body);
    % I don't know why, but the correct rotation matrix seems to be the
    % inverse!!!
    R_body = R_body;
    G_body = [R_body T_body ; zeros(1,3) 1];
    S_body = homo2screw(G_body);

    BodyAng(:,k) = Rot2Joint(R_body);
    
    % My ordering for the twist axes corresponds to the Rxyz Euler angles.
    %
    Theta_Lwing = Rot2Joint(quat2matNEW(q_Lwing_rel));
    Theta_Rwing = Rot2Joint(quat2matNEW(q_Rwing_rel));
    
%     State(:,k) = [T_body
%         q_body
%         q_Lwing_rel
%         q_Rwing_rel];

    State(:,k) = [S_body
        Theta_Lwing
        Theta_Rwing];
end

PAR.params = params;
frames = 1:258;
%% Movie frames
SOLN = State';
for i= 1:length(frames);
    kk = frames(1)+(i-1);

    PAR.imgres = [512 512];
    PAR.imgres_crop = PAR.imgres;
    
    PAR.params.T_Lwing_rel = PAR.params.T_Lwing_rel_All(:,i);
    PAR.params.T_Rwing_rel = PAR.params.T_Rwing_rel_All(:,i);
    
    clear flymod
    [x,y,z] = flymod(SOLN(i,:),PAR.params,PAR);
    
    %[x,y,z] = flymodQ(SOLN(i,:),PAR.params,PAR);
    
    
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    %----------------------
    % Plot Movie Frame or Rendered model image
    %----------------------
    pts = [];
    for j = 1:length(x)
        pts{j} = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
    end
    
       
    for k = 1:PAR.numcam
        
        %Render Model
        [Y,idx,IMout,Nrml] = renderflyMOD(pts,PAR,k);
        file2write = sprintf('%s/%s/cam%03d/%s%06d.bmp', ...
            PAR.videopath,[PAR.filetag 'BW'],k,[PAR.filetag 'BW'],kk);
        imwrite(IMout,file2write,'bmp');
        
    end
end