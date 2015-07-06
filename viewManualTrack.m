% view saved kinematics that have been tracked manually


%clear 
%close all
fignum = 3;
plotscale = .75;

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
%load kine/SavedKinematics/from_gwyneth/exp098_20071209_fixed_eb.mat

% Get Calibration
for j = 1:PAR.numcam
    PAR.DLT(:,j) = data.cal.coeff.(['DLT_' num2str(j)]);
end

load flygenmod
params = flygenmod;

PAR.statedim = 12;

State = zeros(PAR.statedim,data.images.frames);
frames = 1:258;
for k = frames%325:521%data.images.frames
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
    
    %q_Lwing = [1 0 0 0]';
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
    %q_Lwing_rel = [0 0 0 1]';
    
    T_L = [0 0 0];
    %T_L = data.kine.left_wing.data.v_trans(:,k);
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

    %q_Rwing = [1 0 0 0]';
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
    %q_Rwing_rel = [0 0 0 1]';
    T_R = [0 0 0];
    %T_R = data.kine.right_wing.data.v_trans(:,k);
    T_R = T_R(:);

    % Relative translation from body fixed axis to wing joint
    params.T_Rwing_rel_All(:,k) = qxform([-q_body(1:3) ; q_body(4)],T_R - T_body);

%% Calculate the twist representation for these quaternions
    R_body = quat2matNEW(q_body);
    
    G_body = [R_body T_body ; zeros(1,3) 1];
    S_body = homo2screw(G_body);

    BodyAng(:,k) = Rot2Euler(R_body');
    
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
%% Run this code to plot the calculated body angles
%Test the 'predictMotion.m' function
% ff = 94:98;
% [Pnew,tnew] = predictMotion(SOLN(ff,:)');
% ss = [SOLN(ff,:)' Pnew];
% clear BodyAng BodyTrans
% for j = 1:size(ss,2)
%     G = screw2homo(ss(:,j));
%     BodyAng(:,j) = Rot2Joint(G(1:3,1:3));
%     BodyTrans(:,j) = G(1:3,4);
% end
% frames = 1:size(ss,2);
%% Body Coordinates
figure; 
subplot(2,1,1);
%frames = 400:490;%1:size(State,2);
%frames = 1:258;
frames = 325:521;

BodyAng1 = BodyAng;
%BodyAng = data.kine.body.data.eulzyx;
hh = plot(frames,BodyAng(1,frames),'r.-',frames,BodyAng(2,frames),'g.-',frames,BodyAng(3,frames),'b.-');
set(gca,'ytick',-pi:pi/4:pi);
legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
title('Body Angles');
grid

subplot(2,1,2);
hh = plot(frames,BodyTrans(1,frames),'r.-',frames,BodyTrans(2,frames),'g.-',frames,BodyTrans(3,frames),'b.-');
legend(hh,'X','Y','Z');
title('Body Translation');
grid

%% Estimate the process noise values for the linear and angular accleration
% Only use values when fly in taking off
TT = BodyTrans(:,frames);
%AA = State(4:6,170:frames(end));
AA = BodyAng(:,frames);


% Use the 'gradient' function instead of central differencing for the first 
% NN points
NN = 1;
Angveltmp = gradient(AA,PAR.dt);
Angacctmp = gradient(Angveltmp,PAR.dt);
Linveltmp = gradient(TT,PAR.dt);
Linacctmp = gradient(Angveltmp,PAR.dt);

begmid = 1+NN;
endmid = size(AA,2)-NN;

%Calculate accelerations in the middle using central differencing
MidAngAcc = (AA(:,begmid+NN:end) + AA(:,1:endmid-NN) - 2.*AA(:,begmid:endmid))./(NN*PAR.dt)^2;
MidLinAcc = (TT(:,begmid+NN:end) + TT(:,1:endmid-NN) - 2.*TT(:,begmid:endmid))./(NN*PAR.dt)^2;

AngAcc = [Angacctmp(:,1:NN) MidAngAcc Angacctmp(:,end-NN+1:end)];
LinAcc = [Linacctmp(:,1:NN) MidLinAcc Linacctmp(:,end-NN+1:end)];

[mean(AngAcc,2) var(AngAcc,0,2)]
[mean(LinAcc,2) var(LinAcc,0,2)]
%% Wing Coordinates
figure; 
subplot(2,1,1);

hh = plot(frames,State(7,frames),'r.-',frames,State(8,frames),'g.-',frames,State(9,frames),'b.-');
set(gca,'ytick',-pi:pi/4:pi);
legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
title('Left wing');
grid

subplot(2,1,2);
hh = plot(frames,State(10,frames),'r.-',frames,State(11,frames),'g.-',frames,State(12,frames),'b.-');
set(gca,'ytick',-pi:pi/4:pi);
legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
title('Right wing');
grid


%% Plot Joint locations

%params.bodyscale = data.kine.body.data.length(frame) / params.body_tip2tip;
params.bodyscale = mean(data.kine.body.data.length(frames)) / (params.bodylen+params.headlen);
avg_winglength = mean([data.kine.left_wing.data.length(frames)...
    data.kine.right_wing.data.length(frames)]);
params.wingscale = avg_winglength / params.wing_tip2tip;


% Normalize each value by body length
BL = params.bodyscale*(params.bodylen+params.headlen);
RJpts = params.T_Rwing_rel_All(frames) ./ BL;
LJpts = params.T_Lwing_rel_All(frames) ./ BL;

RJTrans = mean(RJpts,2);
LJTrans = mean(LJpts,2);

%Get the mean locations
%make LJ y-coordinate positive to get average absolute coordinate
LJTransNew = LJTrans;
LJTransNew(2) = abs(LJTransNew(2));
mean([RJTrans LJTransNew],2)

rjm = [0.2021 0.1055 -0.1477];
ljm = [0.2021 -0.1055 -0.1477];
RJdelta = RJpts - repmat(rjm',1,size(RJpts,2));
LJdelta = LJpts - repmat(ljm',1,size(LJpts,2));
% JointDisplacement = [LJdelta ; RJdelta];

figure; 
subplot(2,1,1);
hh = plot(frames,LJdelta(1,:),'r*-',frames,LJdelta(2,:),'g*-',frames,LJdelta(3,:),'b*-');
legend(hh,'X','Y','Z');
% quiver3(zeros(1,size(LJdelta,2)),zeros(1,size(LJdelta,2)),zeros(1,size(LJdelta,2)),...
%     LJdelta(1,:),LJdelta(2,:),LJdelta(3,:));
title('Left Joint displacement')

subplot(2,1,2);
hh = plot(frames,RJdelta(1,:),'r.-',frames,RJdelta(2,:),'g.-',frames,RJdelta(3,:),'b.-');
legend(hh,'X','Y','Z');
% quiver3(zeros(1,size(RJdelta,2)),zeros(1,size(RJdelta,2)),zeros(1,size(RJdelta,2)),...
%     RJdelta(1,:),RJdelta(2,:),RJdelta(3,:));
title('Right Joint displacement')

ljmag = sqrt(sum(LJdelta.^2,1));
rjmag = sqrt(sum(RJdelta.^2,1));
figure; 
plot(ljmag,'.-')
hold on
plot(rjmag,'.-')


figure; 
hh = plot(frames,RJpts(1,:),'r.-',frames,RJpts(2,:),'g.-',frames,RJpts(3,:),'b.-');
legend(hh,'X','Y','Z');
hold on; 
plot(frames,RJTrans(1),'r--',frames,RJTrans(2),'g--',frames,RJTrans(3),'b--');

plot(frames,LJpts(1,:),'r*-',frames,LJpts(2,:),'g*-',frames,LJpts(3,:),'b*-');
legend(hh,'X','Y','Z');
hold on; 
plot(frames,LJTrans(1),'r:',frames,LJTrans(2),'g:',frames,LJTrans(3),'b:');
%% Movie frames
SOLN = State';
for i= 1:length(frames);
    kk = frames(1)+(i-1);

    % ------------------------------------
    % Load images
    for cam=1:PAR.numcam
        input_filename = sprintf('%s/%s/cam%03d/%s%06d.bmp', ...
            PAR.videopath,PAR.filetag,cam,PAR.filetag,kk);
        Input(:,:,cam) = imread(input_filename);
    end

    PAR.imgres = size(Input(:,:,1));
    PAR.imgres_crop = PAR.imgres;
    
    PAR.params.T_Lwing_rel = PAR.params.T_Lwing_rel_All(:,i);
    PAR.params.T_Rwing_rel = PAR.params.T_Rwing_rel_All(:,i);
    
    clear flymod
    [x,y,z] = flymod(SOLN(kk,:),PAR.params,PAR);
    
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

    % Initialize projected planar points
    u = cell(PAR.numcam,length(x));
    v = cell(PAR.numcam,length(x));
    uT = cell(PAR.numcam,length(x));
    vT = cell(PAR.numcam,length(x));

    figure(fignum);
    clf;

    im = reshape(Input,PAR.imgres(1),[]);
    imagesc(im);
    colormap gray
    ax1 = gca;
    set(ax1,'units','pixels','visible','off','position',[0 0 PAR.imgres(2)*PAR.numcam PAR.imgres(1)].*plotscale);
    axis image

    %aspectratio = get(ax1,'PlotBoxAspectRatio');
    impos = get(ax1,'position');
    %we want position to be for single image, not composite image
    impos(3) = impos(3)/PAR.numcam;
    color = {'g','r','b'};
    for k = 1:PAR.numcam
        ax2 = axes('position',[0 0 .5 .5]);
        hold on

        DLTparams = PAR.DLT(:,k);
        
        for j = 1:length(pts)
            uv = dlt_3D_to_2D(DLTparams,pts{j});

            u{k,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2));
            v{k,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2));

            uT{k,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2))';
            vT{k,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2))';
            if j == 1
                %Plot Head and tail
                %plot(ax2,uv([1 end],1),uv([1 end],2),'rs');
            end

            hold on;
            patch(u{k,j},v{k,j},[1 1 1],'facecolor','none','edgecolor','w');
            patch(uT{k,j},vT{k,j},[1 1 1],'facecolor','none','edgecolor','w');
        end
        %keyboard
        
        set(ax2,'units','pixels','fontsize',12,'position',impos+[PAR.imgres(2)*(k-1) 0 0 0]*plotscale,'color','none','xlim',...
            [0.5 PAR.imgres(2)+0.5],'ylim',[-0.5 PAR.imgres(1)-0.5],'visible','off',...
            'xdir','normal','ydir','reverse');
        
        %Plot frame # on the 2nd camera image (in the middle)
        if k == 2
            text(PAR.imgres(2)/2 - 50,20,[ '\color{white}Frame ' num2str(kk)], ...
                'fontsize',20);
        end
    end

    figure(fignum);
    set(fignum,'units','pixels','position',[50 150 PAR.imgres(2)*PAR.numcam PAR.imgres(1)].*plotscale);
    %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])
    

    M2(i) = getframe(fignum);
    %M2(kk) = getframe(fignum);
end