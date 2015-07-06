% Initialize the fly model and plot it.
% Changed 11/2/07 to get a better relative wing orientation
clear 
close all

PAR.dt = 1/6000;
PAR.etamax = 0;
PAR.c = 4;
PAR.pixpermm = 1;
PAR.numfly = 1;
PAR.L = 20; %for body
PAR.L2 = 5; %for head
PAR.L3 = 30; %for wing
PAR.T = 15; %theta steps for head and body
PAR.numcam = 3;
PAR.filetag = 'exp098';

load ManualFit_exp098f325
frame = 325;

load flygenmod
params = flygenmod;

%Euler angle rotation functions
Rx = @(theta) [1 0 0
    0 cos(theta) -sin(theta)
    0 sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 -sin(theta)
    0 1 0
    sin(theta) 0  cos(theta)];
Rz = @(theta) [cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];

% Get body params and coordinates
%params.bodyscale = data.kine.body.data.length(frame) / params.body_tip2tip;
params.bodyscale = data.kine.body.data.length(frame) / (params.bodylen+params.headlen);
avg_winglength = mean([data.kine.Left_wing.data.length(frame)
    data.kine.Right_wing.data.length(frame)]);
params.wingscale = avg_winglength / params.wing_tip2tip;

%===============================================================
%Get the body orientation
headpt = data.kine.body.data.coords(1,:,frame);
tailpt = data.kine.body.data.coords(2,:,frame);
coord_diff = headpt - tailpt;
[theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
%theta = data.kine.body.data.angles(1,1,frame);
%phi = data.kine.body.data.angles(1,2,frame);
eul_xyz(1,:) = theta;
eul_xyz(2,:) = phi;
eul_xyz(3,:) = (pi/180)*data.kine.body.data.params(frame);
R_body = Rxyz(eul_xyz);
q_body = dataAngles2quat(eul_xyz);
q_body = q_body(:);


%Calculate inverse rotation to put model in data frame
q_body(2:4) = -q_body(2:4);
R_body = R_body';

% %Get the location of the model head point to determine the tilt that needs
% %to be applied to generative model to match Gwyneth's body frame
% [xx,yy,zz] = flymod([zeros(3,1);zeros(3,1);1;zeros(3,1);1;zeros(3,1);1],params,PAR);
% %I know it lies in the X-Z plane
% model_headpt = [xx{1}(end-1,3) zz{1}(end-1,3)];
% bodytilt = atan2(model_headpt(2),model_headpt(1));
% 
bodytilt = deg2rad(-4);
q_body = quatprodS1(q_body,[cos(bodytilt/2) 0 1*sin(bodytilt/2) 0]);
R_body = R_body*Ry(bodytilt);
%Calculation translation vector -q_body*T10*q_body` + tailpt
%T10 is the approximate location of tailpt in our generative model
%It is the next to last cross-section of model at the dorsal edge
[xbody,ybody,zbody,s,th,X,Frenet,T10] = flybodymod(params.bodyctr,params.bodyrad,params.bodylen,PAR);

%Say T10 is right on the tip of the thorax
T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
T10 = -T10;
[T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),tailpt,q_body);




    
%===============================================================
% Get the wing orientations
% Left
% theta = data.kine.Left_wing.data.angles(1,1,frame);
% phi = data.kine.Left_wing.data.angles(1,2,frame);
base = data.kine.Left_wing.data.coords(1,:,frame);
tip = data.kine.Left_wing.data.coords(2,:,frame);
coord_diff = tip - base;
[theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
eul_xyz(1,:) = theta;
eul_xyz(2,:) = phi;
%eul_xyz(3,:) = (pi/180)*(data.kine.Left_wing.data.params(frame) - 265);
eul_xyz(3,:) = (pi/180)*data.kine.Left_wing.data.params(frame);
q_Lwing = dataAngles2quat(eul_xyz);

R_Lwing = Rxyz(eul_xyz);

%Calculate inverse rotation to put model in data frame
q_Lwing(2:4) = -q_Lwing(2:4);
R_Lwing = R_Lwing';

T_L = base;
T_L = T_L(:);
% I'm going to represent this coordinate transform in the homogeneous
% representation.  
% The homogeneous coordinates are used to calculate the resultant 
% translation, but for some reason, the orientation is off.  Something 
% fishy with quat2mat().  Instead, I will calculate the relative rotations
% for the wings by multiplying the quaternions directly.

RJTrans = [0.52 -0.37 0.18];
LJTrans = [0.52 0.37 0.18];

G_FLprime = [quat2mat(q_Lwing) T_L
    zeros(1,3) 1];
G_FB = [quat2mat(q_body) T_body
    zeros(1,3) 1];
G_BL = [Rx(pi)*Rz(pi/2) LJTrans'
    zeros(1,3) 1];

%This is the local transformation in the wing frame
G_LprimeL = inv(G_FLprime)*G_FB*G_BL;

G_ll = inv(G_FB*G_BL)*G_FLprime;

% in quaternions
atemp = quatprod(q_body,quatprod([cos(pi/2) 1*sin(pi/2) 0 0],[cos(pi/4) 0 0 1*sin(pi/4)]));

%Just take the orientation because I assume that the wing is fixed at the 
%joint. 
%q_Lwing_rel = quat2mat(G_LprimeL(1:3,1:3))';
q_Lwing_rel = quatprod([atemp(1) ; -atemp(2:4)],q_Lwing);
%q_RxRz = quatprod([cos(pi/2) 1*sin(pi/2) 0 0],[cos(pi/4) 0 0 1*sin(pi/4)]);
%q_Lwing_rel = quatMAT([q_Lwing(1) ; -q_Lwing(2:4)])*quatMAT(q_body)*q_RxRz;
T_Lwing_rel = G_LprimeL(1:3,4);

% Right 
% theta = data.kine.Right_wing.data.angles(1,1,frame);
% phi = data.kine.Right_wing.data.angles(1,2,frame);
base = data.kine.Right_wing.data.coords(1,:,frame);
tip = data.kine.Right_wing.data.coords(2,:,frame);
coord_diff = tip - base;
[theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
eul_xyz(1,:) = theta;
eul_xyz(2,:) = phi;
%eul_xyz(3,:) = (pi/180)*(-data.kine.Right_wing.data.params(frame) + 285);
eul_xyz(3,:) = (pi/180)*data.kine.Right_wing.data.params(frame);
q_Rwing = dataAngles2quat(eul_xyz);

%Calculate inverse rotation to put model in data frame
q_Rwing(2:4) = -q_Rwing(2:4);

T_R = base;
T_R = T_R(:);

G_sr = [quat2mat(q_Rwing) T_R
    zeros(1,3) 1];
G_so_body = [quat2mat(q_body) T_body
    zeros(1,3) 1];
G_or_right = [Rz(-pi/2) RJTrans'
    zeros(1,3) 1];

G_sr_ = G_so_body*G_or_right;
G_sr_inv = [G_sr_(1:3,1:3)' -G_sr_(1:3,1:3)'*G_sr_(1:3,4)
    zeros(1,3) 1];

%This is the local transformation in the wing frame
G_rr = G_sr_inv*G_sr;

% in quaternions
atemp = quatprod(q_body,[cos(-pi/4) 0 0 1*sin(-pi/4)]);

%Just take the orientation because I assume that the wing is fixed at the 
%joint. 
%q_Rwing_rel = quat2mat(G_rr(1:3,1:3))';
q_Rwing_rel = quatprod([atemp(1) ; -atemp(2:4)],q_Rwing);
T_Rwing_rel = G_rr(1:3,4);

%===============================================================
%Assemble fly state;
% p = [T_body
%     q_body
%     q_Lwing
%     q_Rwing];

p = [T_body
    q_body
    q_Lwing_rel
    q_Rwing_rel];

%ttemp = [T_L;T_R];
%t_rel = [T_Lwing_rel;T_Rwing_rel];

[x,y,z] = flymodQ(p,params,PAR);
% [x,y,z] = flymod1(p,ttemp,params,PAR);

%=========================================================
%PLOTTING
%=========================================================
% figure(1);
% clf;
% %Plot world reference frame
% ref = eye(3);
% org = zeros(3,3);
% quiver3(org(:,1),org(:,2),org(:,3),ref(:,1),ref(:,2),ref(:,3));
% hold on; 
% 
colors = {'g','r','b'};
for i = 1:length(x)
    surf(x{i},y{i},z{i},'facecolor',colors{i});
    hold on;
end
axis equal
% %view(0,0);

%=========================================================
%PLOTTING on images
fignum = 2;
figure(fignum);
clf;
pts = [];
for i = 1%:length(x)
    pts = [pts ; reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
end
for i = 1:PAR.numcam
    im = imread(['video/' PAR.filetag '/cam00' num2str(i) '/' PAR.filetag ...
        repmat(num2str(0),1,6-(1+floor(log10(frame)))) num2str(frame) '.bmp'],'bmp');
    
    PAR.imgres = size(im);
    
    DLTparams = data.cal.coeff.(['DLT_' num2str(i)]);
    uv = dlt_3D_to_2D(DLTparams,pts);
    
    xax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 11 5 14]);
    yax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 6 14]);
    zax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 5 15]);
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
    
    plot(ax2,uv(:,1),uv(:,2),'.');
    hold on
    %Plot Head and tail
    plot(ax2,uv([1 end],1),uv([1 end],2),'rs');
    
    plot(xax(:,1),xax(:,2),'r.-',yax(:,1),yax(:,2),'g.-',zax(:,1),zax(:,2),'b.-')
    %label the axes
    text(xax(2,1),xax(2,2),'X','fontsize',12)
    text(yax(2,1),yax(2,2),'Y','fontsize',12)
    text(zax(2,1),zax(2,2),'Z','fontsize',12)
 
    set(ax2,'fontsize',12,'position',impos,'color','none','xlim',...
        [0.5 PAR.imgres(2)+0.5],'ylim',[-0.5 PAR.imgres(1)-0.5],'visible','on',...
        'xdir','normal','ydir','reverse');
    set(fignum,'position',[0 500 1280 427])
        
end
    
    

