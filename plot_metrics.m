%plot_metrics.m
%
%plot_metrics plots the state variables of a specified sequence over a
%specified number of frames.  
%It will plot the estimated values and/or values from a prior database.
%
%clear all
%close all

PAR.videopath = '../FlyVideo/';
PAR.filetag = 'exp098000';

PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.stub = [PAR.filetag];
fignum = 1;

% --------------------------------
% Load the ManualFit Parameters
load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);

%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
frames = [325 521]; % number of frames in movie

samplerate = 1;
movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);

interptag = 1;
%--------------------------------
% Define the Tracking Parameters
PAR.pixpermm = 1;
PAR.numfly = 1;
%Number of parameters of the model (i.e. 8 control points)
PAR.mdlpar = 15*ones(1,PAR.numfly);
PAR.statedim = PAR.mdlpar;
PAR.modelfun_H = @modcurvesplineP;
PAR.etamax = 0;

%spline order
PAR.c = 4;
PAR.L1 = 15; %# of steps for body along length
PAR.L2 = 6; %# of steps for head along length
PAR.L3 = 25; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;
solidx = [1 frames(2)];

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;


soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
SOLN = zeros(length(solidx(1):solidx(2)),PAR.numfly*PAR.statedim);
if interptag == 1
    for i=1:length(movidx)
        load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(movidx(i)) ...
            '.mat']);
%         load([PAR.solutionpath 'fly_' PAR.stub '_BodyPredictQ/fly' num2str(movidx(i)) '.mat']);
        
        soln1(i,:) = xh(1:PAR.numfly*PAR.statedim)';
        clear xh InternalVariablesDS
    end
end

%SOLN = soln1;

%Now, perform interpolation to calculate the state at in between frames
for k = 1:size(soln1,2)
    SOLN(frames(1):frames(2),k) = interp1(movidx,soln1(:,k),frames(1):frames(2),'spline')';
end

%% ---------------------------------------------------------
% Now, Calculate the spatial velocity.
% I'm interested in the angular components only to model the 
% Process noise corectly in the motion model
% scr = SOLN(frames(1):frames(2),1:6);
% 
% %=====================================
% % calculate theta-dot
% % For twist representation only
% %======================================
% theta_vec = sqrt(sum(scr(:,4:6).^2,2));
% theta_dot = gradient(theta_vec,PAR.dt);
% 
% % Use MSE quintic spline method from Walker, JEB 1998 to calculate
% % displacements, velocities, and accelerations.  It uses the Matlab spline 
% % toolbox 
% 
% time = (0:size(scr,1)-1)*PAR.dt*1000;
% 
% % The expected error for correctly locating an pixel location is
% tol = sqrt(2*0.5^2);  %see Walker equation #13 with h = 0;
% % Convert to body length units
% len = 2.2631;
% pixpermm = 38.8754;
% tol = tol*(1/pixpermm)*(1/sum(len));
% % data point weighting
% w = ones(size(time));
% 
% timeSec = time/1000;
% timeSec = timeSec(:);%make column vector
% 
% % =========================================
% % Calculations for the splined screw magnitude, theta
% [thetafun,theta_vecNew] = spaps(timeSec,theta_vec,tol,w,3);
% thetadotfun = fnder(thetafun,1);
% 
% %Calculate Velocity 
% theta_dotNew = fnval(thetadotfun,timeSec)';
% 
% %Use the splined values instead
% theta_vec = theta_vecNew;
% theta_vec = theta_vec(:);%make column vector
% 
% theta_dot = theta_dotNew;
% theta_dot = theta_dot(:);%make column vector
% 
% 
% skew = @(p) [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];
% 
% %Make the angular part unit length
% scr = scr ./ [ones(size(scr,1),3) repmat(theta_vec,1,3)];
% 
% %Scale by theta dot to get spatial velocity
% scr = scr.* repmat(theta_dot,1,6);
% 
% ang_vel = scr(:,4:6);
% for i = 1:3
%     ang_acc(:,i) = gradient(ang_vel(:,i),PAR.dt);
% end
% ang_var = var(ang_acc,0,1);
% 
% %---------------------------------------------------------
% % Now, Calculate the body velocity.
% % I'm interested in the linear components (velocity of body frame
% % relative to spatial frame) only to model the 
% % Process noise corectly in the motion model
% scr = SOLN(frames(1):frames(2),1:6);
% 
% % Take the inverse screw motion for the body velocity
% scr = -scr;
% 
% %Make the angular part unit length
% scr = scr ./ [ones(size(scr,1),3) repmat(theta_vec,1,3)];
% 
% %Scale by theta dot to get body velocity
% scr = scr.* repmat(theta_dot,1,6);
% 
% lin_vel = scr(:,1:3);
% for i = 1:3
%     lin_acc(:,i) = gradient(lin_vel(:,i),PAR.dt);
% end
% lin_var = var(lin_acc,0,1);


%% ===============================================================

% Convert the automatically captured body data into Euler components
BodyTrans_auto = zeros(3,frames(2));
BodyAng_auto = zeros(3,frames(2));
for k = frames(1):frames(2)
    %G = screw2homo(SOLN(k,1:6));
    Rbod = quat2matNEW(SOLN(k,4:7));
    BodyAng_auto(:,k) = Rot2Euler(Rbod');%Rot2Euler(G(1:3,1:3)');
    BodyTrans_auto(:,k) = SOLN(k,1:3)';%G(1:3,4);
end

SOLN = SOLN';

%% Load the Manually Tracked data
%load kine/SavedKinematics/Wings/exp083_forEbraheem_corAlpha_splined_adj.mat
load kine/SavedKinematics/Wings/exp101_ebraheem.mat
%load kine/SavedKinematics/from_gwyneth/exp098_20071209_fixed_eb.mat
%load([PAR.solutionpath 'fly_' PAR.stub '/KineManualFit_exp098f325.mat'])
frames1 = [1 2];
%State = zeros(PAR.statedim,length(frames1));

BodyTrans = zeros(3,frames1(2));
BodyAng = zeros(3,frames1(2));

for k = frames1
    
%% Get the body orientation
    q_body = data.kine.body.data.quat(:,k);
    
    %I store the scalar part of quaternion at the end
    q_body = [q_body(2:4); q_body(1)];
    
    % add a rotation by pi along the roll axis since our model has the body
    % fixed frame with z-axis pointing ventral.
    q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);
    
    % Get Translation
    T_body = data.kine.body.data.v_trans(:,k);
    %Calculation translation vector -q_body*T10*q_body` + T_body
    %T10 is the approximate location of tailpt in our generative model
    %It is the next to last cross-section of model at the dorsal edge
    [xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(ManualFit.params.bodyctr,...
        ManualFit.params.bodyrad,ManualFit.params.bodylen,PAR);

    %Say T10 is right on the tip of the thorax
    %T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
    %T10 = -T10;
    %[T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),T_body,q_body);
    
    %Calculate the equivalent location using manually fit body length as the 
    %T_body used in our model for proper error comparison.
    bodyvec = diff(data.kine.body.data.coords(:,:,k),1,1);
    bodyvec = -bodyvec./norm(bodyvec);
    T_body = T_body + .4226*data.kine.body.data.length(k).*bodyvec';    
    

    BodyTrans(:,k) = T_body;
    
%% Get the Left wing orientation
    %q_Lwing = [0 0 0 1]';
    q_Lwing = data.kine.left_wing.data.quat(:,k);
    
    %I store the scalar part of quaternion at the end
    q_Lwing = [q_Lwing(2:4); q_Lwing(1)];

    % Premultiply rotation by this alignment quaternion that takes into account
    % the orientation of coordinate axis fixed to the left wing.
    q_Lwingaxisalign = quat2matNEW([0 -1 0;1 0 0;0 0 1]);
    q_Lwing = quatprod(q_Lwing,q_Lwingaxisalign);

    % I will calculate the relative rotations
    % for the wings by multiplying the quaternions directly.
    % Just take the orientation because I assume that the wing is fixed at the
    % joint.
    q_Lwing_rel = quatprod([-q_body(1:3) ; q_body(4)],q_Lwing);
    
    T_L = [0 0 0];
    %T_L = data.kine.left_wing.data.v_trans(:,k);
    T_L = T_L(:);
   
    % Relative translation from body fixed axis to wing joint
    params.T_Lwing_rel_All(:,k) = qxform([-q_body(1:3) ; q_body(4)],T_L - T_body);

%% Right Wing
    %q_Rwing = [0 0 0 1]';
    q_Rwing = data.kine.right_wing.data.quat(:,k);
    %I store the scalar part of quaternion at the end
    q_Rwing = [q_Rwing(2:4); q_Rwing(1)];

    % Premultiply rotation by this alignment quaternion that takes into account
    % the orientation of coordinate axis fixed to the right wing.
    q_Rwingaxisalign = quat2matNEW([0 1 0;1 0 0;0 0 -1]);
    q_Rwing = quatprod(q_Rwing,q_Rwingaxisalign);

    %Just take the orientation because I assume that the wing is fixed at the
    %joint.
    q_Rwing_rel = quatprod([-q_body(1:3); q_body(4)],q_Rwing);
    
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
    
    StateQ(:,k) = [T_body
        q_body
        q_Lwing_rel
        q_Rwing_rel];

    State(:,k) = [S_body
        Theta_Lwing
        Theta_Rwing];
end

%Common frames between manual and automatically tracked data
%ff = 90:425;
ff = frames(1):frames(2);
%% Body Coordinates
figure; 
subplot(1,2,1);
mkr = {'--','-'};

BodyAng = BodyAng.*(180/pi);
BodyAng_auto = BodyAng_auto.*(180/pi);
tt = (0:length(ff)-1)*PAR.dt*1000;

hh = [];

%plot prior database
%hh = plot(tt,BodyAng(1,ff),['r' mkr{1}],tt,BodyAng(2,ff),['g' mkr{1}],tt,BodyAng(3,ff),['b' mkr{1}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',1)
end

hold on;
%plot tracked data
hh = plot(tt,BodyAng_auto(1,ff),['r' mkr{2}],tt,BodyAng_auto(2,ff),['g' mkr{2}],tt,BodyAng_auto(3,ff),['b' mkr{2}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',2)
end

%set(gca,'ytick',-pi:pi/4:pi,'xlim',[ff(1) ff(end)]);
set(gca,'xlim',[tt(1) tt(end)]);
legend(hh,'roll','pitch','yaw');
title('Body Angles');
ylabel('degrees');
xlabel('ms');
grid

subplot(1,2,2);
%plot prior database
%hh = plot(tt,BodyTrans(1,ff),['r' mkr{1}],tt,BodyTrans(2,ff),['g' mkr{1}],tt,BodyTrans(3,ff),['b' mkr{1}]);
for i = 1:length(hh)
    set(hh(i),'linewidth',1)
end

hold on;
%plot tracked data
hh = plot(tt,BodyTrans_auto(1,ff),['r' mkr{2}],tt,BodyTrans_auto(2,ff),['g' mkr{2}],tt,BodyTrans_auto(3,ff),['b' mkr{2}]);
for i = 1:length(hh)
    set(hh(i),'linewidth',2)
end

set(gca,'xlim',[tt(1) tt(end)])
legend(hh,'X','Y','Z');
title('Body Translation');
ylabel('mm');
xlabel('ms');
grid


%% Wing Coordinates
figure; 
% --------------------------------
% Left wing
subplot(2,1,1);
%plot prior database
%hh = plot(ff,State(7,ff),['r' mkr{1}],ff,State(8,ff),['g' mkr{1}],ff,State(9,ff),['b' mkr{1}]);

hold on;
%plot tracked data
%  hh = plot(ff,SOLN(7,ff),['r' mkr{2}],ff,SOLN(8,ff),['g' mkr{2}],ff,SOLN(9,ff),['b' mkr{2}]);
hh = plot(ff,SOLN(8,ff),['r' mkr{2}],ff,SOLN(9,ff),['g' mkr{2}],...
    ff,SOLN(10,ff),['b' mkr{2}],ff,SOLN(11,ff),['k' mkr{2}]);

%set(gca,'ytick',-pi:pi/4:pi);
%legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
legend(hh,'q_1','q_2','q_3','q_4');
title('Left wing');
grid

% --------------------------------
% Right wing
subplot(2,1,2);
%plot prior database
%hh = plot(ff,State(10,ff),['r' mkr{1}],ff,State(11,ff),['g' mkr{1}],ff,State(12,ff),['b' mkr{1}]);

hold on;
%plot tracked data
%hh = plot(ff,SOLN(10,ff),['r' mkr{2}],ff,SOLN(11,ff),['g' mkr{2}],ff,SOLN(12,ff),['b' mkr{2}]);
hh = plot(ff,SOLN(12,ff),['r' mkr{2}],ff,SOLN(13,ff),['g' mkr{2}],...
    ff,SOLN(14,ff),['b' mkr{2}],ff,SOLN(15,ff),['k' mkr{2}]);

% set(gca,'ytick',-pi:pi/4:pi);
% legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
legend(hh,'q_1','q_2','q_3','q_4');
title('Right wing');
grid