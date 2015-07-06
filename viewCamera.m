clear all
close all

load kine/Cali'brations for Kine'/cal_for_20031217_done_20040205.mat

PAR.dt = 1/6000;
PAR.etamax = 0;
PAR.c = 4;
PAR.pixpermm = 1;
PAR.numfly = 1;
PAR.L1 = 50; %# of steps for body along length
PAR.L2 = 20; %# of steps for head along length
PAR.L3 = 30; %# of steps for wing around the boundary
PAR.T1 = 25; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

PAR.numcam = 3;
PAR.filetag = 'exp083';
PAR.imgres = [512 512];

% Convert camera coordinates so that origin is at the lower left with
% z-axis pointing towards you
% L_1(:,2) = PAR.imgres(2) - L_1(:,2);
% L_2(:,2) = PAR.imgres(2) - L_2(:,2);
% L_3(:,2) = PAR.imgres(2) - L_3(:,2);

L(:,:,1) = L_1;
L(:,:,2) = L_2;
L(:,:,3) = L_3;

for i = 1:3
    a(:,i) = dltfu_iter(F,L(:,:,i));
    %a(:,i) = dltfu_nonlin(F,L(:,:,i));
    cam(i) = dlt2cam(a(:,i));
end


load ManualFit_exp083f415
%load kine/SavedKinematics/Wings/exp101_ebraheem.mat

frame = 415;

load flygenmod
params = flygenmod;

%Euler angle rotation functions
Rx = @(theta) [1 0 0
    0 cos(theta) -sin(theta)
    0 sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 sin(theta)
    0 1 0
    -sin(theta) 0  cos(theta)];
Rz = @(theta) [cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];

%===============================================================
%% Get body params and coordinates
%params.bodyscale = data.kine.body.data.length(frame) / params.body_tip2tip;
params.bodyscale = data.kine.body.data.length(frame) / (params.bodylen+params.headlen);
avg_winglength = mean([data.kine.left_wing.data.length(frame)
    data.kine.right_wing.data.length(frame)]);
params.wingscale = avg_winglength / params.wing_tip2tip;

%===============================================================
%% Get the body orientation
q_body = data.kine.body.data.quat(:,frame);

%I store the scalar part of quaternion at the end
q_body = [q_body(2:4); q_body(1)];

% add a rotation by pi along the roll axis since our model has the body
% fixed frame with z-axis pointing ventral.
q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);

T_body = data.kine.body.data.v_trans(:,frame);

%Calculation translation vector -q_body*T10*q_body` + T_body
%T10 is the approximate location of tailpt in our generative model
%It is the next to last cross-section of model at the dorsal edge
[xbody,ybody,zbody,s,th,X,Frenet,T10] = flybodymod(params.bodyctr,params.bodyrad,params.bodylen,PAR);

%Say T10 is right on the tip of the thorax
T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
T10 = -T10;
[T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),T_body,q_body);

    
%===============================================================
% Get the wing orientations
%% Left
q_Lwing = data.kine.left_wing.data.quat(:,frame);

%I store the scalar part of quaternion at the end
q_Lwing = [q_Lwing(2:4); q_Lwing(1)];

% Premultiply rotation by this alignment quaternion that takes into account
% the orientation of coordinate axis fixed to the left wing.
q_Lwingaxisalign = quat2matNEW([0 -1 0;1 0 0;0 0 1]);
q_Lwing = quatprod(q_Lwing,q_Lwingaxisalign);


T_L = data.kine.left_wing.data.v_trans(:,frame);
T_L = T_L(:);


% I will calculate the relative rotations
% for the wings by multiplying the quaternions directly.
% Just take the orientation because I assume that the wing is fixed at the 
% joint. 
q_Lwing_rel = quatprod([-q_body(1:3) ; q_body(4)],q_Lwing);

% Relative translation from body fixed axis to wing joint
params.T_Lwing_rel = qxform([-q_body(1:3) ; q_body(4)],T_L - T_body);


%% Right Wing
q_Rwing = data.kine.right_wing.data.quat(:,frame);

%I store the scalar part of quaternion at the end
q_Rwing = [q_Rwing(2:4); q_Rwing(1)];

% Premultiply rotation by this alignment quaternion that takes into account
% the orientation of coordinate axis fixed to the right wing.
q_Rwingaxisalign = quat2matNEW([0 1 0;1 0 0;0 0 -1]);
q_Rwing = quatprod(q_Rwing,q_Rwingaxisalign);


T_R = data.kine.right_wing.data.v_trans(:,frame);
T_R = T_R(:);

%Just take the orientation because I assume that the wing is fixed at the 
%joint. 
q_Rwing_rel = quatprod([-q_body(1:3); q_body(4)],q_Rwing);

% Relative translation from body fixed axis to wing joint
params.T_Rwing_rel = qxform([-q_body(1:3) ; q_body(4)],T_R - T_body);
%===============================================================
%% Assemble fly state;
% p = [T_body
%     q_body
%     q_Lwing
%     q_Rwing];

pQ = [T_body
    q_body
    q_Lwing_rel
    q_Rwing_rel];

%ttemp = [T_L;T_R];
%t_rel = [T_Lwing_rel;T_Rwing_rel];

[x,y,z] = flymodQ(pQ,params,PAR);
%[x,y,z] = flymodQ([zeros(1,3) zeros(1,3) 1 zeros(1,3) 1 zeros(1,3) 1],params,PAR);
% [x,y,z] = flymod1(p,ttemp,params,PAR);


%% Calculate the twist representation for these quaternions
R_body = quat2matNEW(q_body);

G_body = [R_body T_body ; zeros(1,3) 1];
S_body = homo2screw(G_body);

% My ordering for the twist axes corresponds to the Rxyz Euler angles.
% 
Theta_Lwing = Rot2Joint(quat2matNEW(q_Lwing_rel));
Theta_Rwing = Rot2Joint(quat2matNEW(q_Rwing_rel));

p = [S_body
    Theta_Lwing
    Theta_Rwing];

clear flymod
[x,y,z] = flymod(p,params,PAR);
for j = 1:length(x);
    PAR.modsample(j) = size(x{j},1);
end
%=========================================================

%=========================================================
%% PLOTTING
%=========================================================
% figure(3);
% clf;

%Find optimal viewing angle and point camera is looking at
cam = FindOptViewAngle(a,L,F,cam,PAR);
    
colors = {'g','r','b'};
flyrenIM = zeros([PAR.imgres PAR.numcam]);
for i = 1:PAR.numcam
    %subplot(1,3,i);
    figure(i);
    clf;
    
    %Plot world reference frame
    ref = eye(3);
    org = zeros(3,3);
    %quiver3(org(:,1),org(:,2),org(:,3),ref(:,1),ref(:,2),ref(:,3));
    %hold on;
    
    for n = 1:length(x)
        surf(x{n},y{n},z{n},'facecolor',colors{n});
        hold on;
    end

    axis equal
    ax(i) = gca;
    
    set(ax(i),'CameraPosition',cam(i).C,'CameraUpVector',cam(i).R(2,:),...
        'Projection','perspective','CameraTarget',cam(i).xyz,'CameraViewAngle',...
        cam(i).viewangle);
    
    %flyrenIM(:,:,i) = renderfly(i,PAR.imgres);
    
end

Flypts = cell(3,1);
R = [];
for n = 1:length(x)
    Flypts{n} = [reshape(x{n},[],1) reshape(y{n},[],1) reshape(z{n},[],1)];
    C(:,n) = cam(n).C;
    R = [R ; cam(n).R];
end
Flypts = cell2mat(Flypts)';
Flypts = [Flypts ; zeros(1,size(Flypts,2))];

drawscene(Flypts,C,R,4,'cloud','Fly Experimental Setup');
% %view(0,0);

%=========================================================
%PLOTTING on images
fignum = 5;
figure(fignum);
clf;
pts = [];
for i = 1:length(x)
    pts{i} = [reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
end

% Initialize projected planar points
u = cell(PAR.numcam,length(x));
v = cell(PAR.numcam,length(x));
uT = cell(PAR.numcam,length(x));
vT = cell(PAR.numcam,length(x));

for i = 1:PAR.numcam
    im = imread(['video/' PAR.filetag '/cam00' num2str(i) '/' PAR.filetag ...
        repmat(num2str(0),1,6-(1+floor(log10(frame)))) num2str(frame) '.bmp'],'bmp');
    
    PAR.imgres = size(im);
    
    DLTparams = data.cal.coeff.(['DLT_' num2str(i)]);
    DLT(:,i) = DLTparams;
    
    xax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 11 5 14]);
    yax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 6 14]);
    zax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 5 15]);
    %figure(10+i);
    subplot(1,3,i);
    imagesc(im);
    colormap gray

    ax1 = gca;
    set(ax1,'visible','off');
    axis tight
    
    aspectratio = get(ax1,'PlotBoxAspectRatio');
    impos = get(ax1,'position');
    ax2 = axes('position',[0 0 .5 .5],'PlotBoxAspectRatio',aspectratio);
    
    subax(2*i-1) = ax1;
    subax(2*i) = ax2;
    
    for j = 1:length(pts)
        uv = dlt_3D_to_2D(DLTparams,pts{j});
        
        u{i,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2));
        v{i,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2));
        
        uT{i,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2))';
        vT{i,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2))';
        
        if j == 1
            %Plot Head and tail
            plot(ax2,uv([1 end],1),uv([1 end],2),'rs');
        end
    
        hold on;
        patch(u{i,j},v{i,j},[1 1 1],'facecolor','none','edgecolor','w');
        patch(uT{i,j},vT{i,j},[1 1 1],'facecolor','none','edgecolor','w');
    end
    
    plot(xax(:,1),xax(:,2),'r.-',yax(:,1),yax(:,2),'g.-',zax(:,1),zax(:,2),'b.-')
    %label the axes
    text(xax(2,1),xax(2,2),'X','fontsize',12)
    text(yax(2,1),yax(2,2),'Y','fontsize',12)
    text(zax(2,1),zax(2,2),'Z','fontsize',12)
 
    set(ax2,'fontsize',12,'position',impos,'color','none','xlim',...
        [0.5 PAR.imgres(2)+0.5],'ylim',[-0.5 PAR.imgres(1)-0.5],'visible','on',...
        'xdir','normal','ydir','reverse');
    set(fignum,'position',[0 500 3*512*.75 512*.75])
        
end

% ManualFit.frame = frame;
% ManualFit.soln = p;
% ManualFit.solnQ = pQ;
% ManualFit.DLT = DLT;
% ManualFit.params = params;
% ManualFit.cam = cam;
% save(['ManualFit_' PAR.filetag],'ManualFit');