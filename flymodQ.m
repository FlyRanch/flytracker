function [x,y,z] = flymodQ(p,params,PAR)
% function [x,y,z] = flymod(p,params,PAR)
%
% x,y,z are [1 3] cells that contain the body/head points, left wing, and
% right wing points, respectively;
%
% p is the fly state denoted as
% p = [ T_body      - global translation to center of fly body
%       q_body      - quaternion for body orientation
%       q_Lwing     - quaternion for left wing orientation
%       q_Rwing];   - quaternion for right wing orientation
%
% params is a structure that has all of the model parameters that encode
% the shape.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!
% T_body goes to center of fly body (i.e. thorax), not center of body and
% head model assembled
% % bodylen + headlen = 2.3276 + 0.4265 = 2.7541
% 0.5*bodylen / 2.7541 = .4226
% So T_body is located at 42.25% along the fly's body axis from the tail.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PAR.dt = 1/6000;
% PAR.etamax = 0;
% PAR.c = 4;
% PAR.pixpermm = 1;
% PAR.numfly = 1;
% PAR.L = 40;

persistent xorg yorg zorg

x = cell(3,1);
y = x;
z = x;

% the posture transformations
Tbody = p(1:3);
qbody = p(4:7);

BodyOnly = length(p) < 8;
if ~BodyOnly
    qLwing = p(8:11);
    qRwing = p(12:15);
end

% Currently, disable the fly joint location estimation
if length(p) < 0
    DeltaLJ = p(16:18);
    DeltaRJ = p(19:21);
else
    DeltaLJ = [0 0 0]';
    DeltaRJ = [0 0 0]';
end

if isempty(xorg)
    %If these points don't exist from a previous function call, create them

    %=========================================================
    %BODY
    %=========================================================
    %Instead of using the bent fly model from W. Dickson, let's just assume
    %that the body axis is a straight line. Look inside 'flybodymod' for
    %manipulation.
    [xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(params.bodyctr,params.bodyrad,params.bodylen,PAR);

    %=========================================================
    %HEAD
    %=========================================================
    %Transformation to the neck;
    % headtilt = 18*(pi/180);
    headtilt = 0;
    Ry = [cos(headtilt) 0 -sin(headtilt)
        0 1 0
        sin(headtilt) 0 cos(headtilt)];
    %Rneck = Ry*[Frenet.T(end,:)' Frenet.N(end,:)' Frenet.B(end,:)'];
    Rneck = Ry*[Frenet.T(end,:)' -Frenet.B(end,:)' Frenet.N(end,:)'];

    %Tneck = X(end,:) + [-0.020 0 0]; %Shift head back a little
    Tneck = X(end,:);
   %Transformation to put reference frame at bottom of head and shift
    % Thead = [params.headlen/2 0 -0.1];
    Thead = [params.headlen/2 0 0];

    [xhead,yhead,zhead] = flyheadmod(params.headrad,params.headlen,PAR);

    %Translate reference frame at bottom of head
    [xhead,yhead,zhead] = xform_surf(xhead,yhead,zhead,Thead);

    %Transform to neck of body
    [xhead,yhead,zhead] = xform_surf(xhead,yhead,zhead,Tneck,Rneck);

    %Put head and body points together.
    x{1} = [xbody ; xhead];
    y{1} = [ybody ; yhead];
    z{1} = [zbody ; zhead];

    %=========================================================
    %WING
    %=========================================================
    %[x{2},y{2},z{2}] = flywingmod(params.wing,params.winglen,PAR);
    [x{2},y{2},z{2}] = flywingmodL(params.wing,params.winglen,PAR);
    
    %keyboard
    
    %Make Right Wing too
    % x{3} = x{2};
    % y{3} = y{2};
    % z{3} = z{2};
    [x{3},y{3},z{3}] = flywingmodR(params.wing,params.winglen,PAR);

    xorg = x;
    yorg = y;
    zorg = z;
else
    % these points already exist, so just assign them to the appropriate
    % variables
    x = xorg;
    y = yorg;
    z = zorg;
end

if ~BodyOnly
    %Approximate translation from center of fly body to the wing joint location
    % JTrans = [0 +/-0.35 0.5]

    %According to Gwyneth's fly frame model, the fly's joints are approximately
    %2/3 along length from tail to head (This coresponds to 1/6 from midpoint
    %towards head). 20% of length out from the centerline
    %and 16% up from centerline towards dorsal side.

    % RJTrans = params.bodyscale*((params.bodylen+params.headlen)*[1/6 0.1 -0.1]...
    %     + [params.headlen/2 0 0]);
    % LJTrans = params.bodyscale*((params.bodylen+params.headlen)*[1/6 -0.1 -0.1]...
    %     + [params.headlen/2 0 0]);

    % These values are calculated from the hand tracked data 'exp101_ebraheem'
    % See 'viewManualTrack' for details of calculation
    BL = params.bodyscale*(params.bodylen+params.headlen);
    RJTrans = BL.*([0.2021 0.1055 -0.1477] + DeltaRJ');
    LJTrans = BL.*([0.2021 -0.1055 -0.1477] + DeltaLJ');
    % RJTrans = [0 0 0];
    % LJTrans = [0 0 0];
    % LJTrans = params.T_Lwing_rel;
    % RJTrans = params.T_Rwing_rel;

    %Rotate about z-axis by 180 degrees
    Rxpi = [1 0 0
        0 -1 0
        0 0 -1];

    Rzpi2 = [0 1 0
        -1 0 0
        0 0 1];

    Rzmpi2 = [0 1 0
        -1 0 0
        0 0 1];

    % T_L = [0 0 0];
    % T_R = [0 0 0];
    %T_L = params.T_Lwing_rel;
    %T_R = params.T_Rwing_rel;
    % T_L = tt(1:3);
    % T_R = tt(4:6);
    %keyboard
    % Now apply wing rotations about the joint
    %xformq_surfMAT(x{2},y{2},z{2},LJTrans,qLwing,params.wingscale);
    [x{2},y{2},z{2}] = xformq_surf(x{2},y{2},z{2},LJTrans,qLwing,params.wingscale);
    [x{3},y{3},z{3}] = xformq_surf(x{3},y{3},z{3},RJTrans,qRwing,params.wingscale);

    % Transform Left wing to frame 'F'
    %[x{2},y{2},z{2}] = xform_surf(x{2},y{2},z{2},LJTrans,Rxpi*Rzpi2);
    %[x{2},y{2},z{2}] = xform_surf(x{2},y{2},z{2},LJTrans);

    % Transform Right wing to frame 'F'
    %[x{3},y{3},z{3}] = xform_surf(x{3},y{3},z{3},RJTrans,Rzmpi2);
    %[x{3},y{3},z{3}] = xform_surf(x{3},y{3},z{3},RJTrans);
end

%=========================================================
%=========================================================

if BodyOnly
    x = x(1);
    y = y(1);
    z = z(1);
end

% global transformation for all body parts (and scaling!)
for k = 1:length(x)
    if k == 1
        [x{k},y{k},z{k}] = xformq_surf(x{k},y{k},z{k},Tbody,qbody,params.bodyscale);
    else
        %keyboard
        [x{k},y{k},z{k}] = xformq_surf(x{k},y{k},z{k},Tbody,qbody,1);
    end
end


