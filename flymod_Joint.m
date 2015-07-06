function [x,y,z] = flymod_Joint(p,params,PAR)
% function [x,y,z] = flymod_Joint(p,params,PAR)
%
% x,y,z are [3 1] cells that contain the body/head points, left wing, and
% right wing points, respectively;
% 
% Twist notation is taken from Murray Book "Robotic Manipulation"
%
% p is the fly state denoted as
% p = [ xi          - global screw motion of the fly model xi = [v ; w]
%       theta_Lwing     - [t1 t2 t3] for left wing orientation, the
%                          magnitude of rotation in radians around each of 
%                          the wing's rotational axes
%       theta_Rwing     - [t1 t2 t3] for left wing orientation, the
%                          magnitude of rotation in radians around each of 
%                          the wing's rotational axes
%       DeltaLJ         - Left joint displacement from nominal value [3 1]
%                         vector
%       DeltaRJ         - Right joint displacement from nominal value [3 1]
%                         vector
% params is a structure that has all of the model parameters that encode
% the shape.

persistent xorg yorg zorg

x = cell(3,1);
y = x;
z = x;

p = p(:);
% the posture transformations
bodyscrew = p(1:6);
Gbody = screw2homo(bodyscrew);

Tbody = Gbody(1:3,4);
Rbody = Gbody(1:3,1:3);

% joint angles for the wings.
% Theta(1:3) ==> Left wing
% Theta(4:6) ==> Left wing
Theta = p(7:12);

DeltaLJ = p(13:15);
DeltaRJ = p(16:18);
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
    % headtilt = -18*(pi/180);
    headtilt = 0;
    Ry = [cos(headtilt) 0 sin(headtilt)
        0 1 0
        -sin(headtilt) 0 cos(headtilt)];
    Rneck = Ry*[Frenet.T(end,:)' Frenet.N(end,:)' Frenet.B(end,:)'];

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

%Approximate translation from center of fly body to the wing joint location
% JTrans = [0 +/-0.35 0.5]

%According to Gwyneth's fly frame model, the fly's joints are approximately
%2/3 along length from tail to head (This coresponds to 1/6 from midpoint 
%towards head). 20% of length out from the centerline
%and 16% up from centerline towards dorsal side.

% RJTrans = params.bodyscale*((params.bodylen+params.headlen)*[1/6 0.1 -0.15]...
%     + [params.headlen/2 0 0]);
% LJTrans = params.bodyscale*((params.bodylen+params.headlen)*[1/6 -0.1 -0.15]...
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

%Rotate about x-axis by 180 degrees
Rxpi = [1 0 0
    0 -1 0
    0 0 -1];

%Rotation about z by pi/2 (counter-clockwise)
Rzpi2 = [0 1 0
    -1 0 0
    0 0 1];

%Rotation about z by -pi/2 (clockwise)
Rzmpi2 = [0 1 0
    -1 0 0
    0 0 1];

T_L = [0 0 0];
T_R = [0 0 0];
% T_L = tt(1:3);
% T_R = tt(4:6);

% Transform Left wing to frame 'L'
%[x{2},y{2},z{2}] = xform_surf(x{2},y{2},z{2},LJTrans);

% Transform Right wing to frame 'R'
%[x{3},y{3},z{3}] = xform_surf(x{3},y{3},z{3},RJTrans);

% Now apply wing rotations about the joint according to joint angles
%
% Twist vectors for the Left wing
% S(:,1) = [0 LJTrans(3) -LJTrans(2) 1 0 0]'; 
% S(:,2) = [-LJTrans(3) 0 LJTrans(1) 0 1 0]';
% S(:,3) = [LJTrans(2) -LJTrans(1) 0 0 0 1]';

% and Right wing
% S(:,4) = [0 RJTrans(3) -RJTrans(2) 1 0 0]'; 
% S(:,5) = [-RJTrans(3) 0 RJTrans(1) 0 1 0]';
% S(:,6) = [RJTrans(2) -RJTrans(1) 0 0 0 1]';

S = [zeros(3,6) ; [eye(3) eye(3)]];
%Scale the twists according to joint angle
S(4:6,:) = S(4:6,:).*repmat(Theta',3,1);

G_BL = screw2homo(S(:,1))*screw2homo(S(:,2))*screw2homo(S(:,3));
G_BR = screw2homo(S(:,4))*screw2homo(S(:,5))*screw2homo(S(:,6));

%The translation part of the homogeneous transformations above will always
%be zero, so just include the wing joint translations below
%keyboard
[x{2},y{2},z{2}] = xform_surf(x{2},y{2},z{2},LJTrans,G_BL(1:3,1:3),params.wingscale);
[x{3},y{3},z{3}] = xform_surf(x{3},y{3},z{3},RJTrans,G_BR(1:3,1:3),params.wingscale);


%=========================================================
%=========================================================

% global transformation for all body parts (and scaling!)
for k = 1:length(x)
    if k == 1 %body
        [x{k},y{k},z{k}] = xform_surf(x{k},y{k},z{k},Tbody,Rbody,params.bodyscale);
    else %wings
        [x{k},y{k},z{k}] = xform_surf(x{k},y{k},z{k},Tbody,Rbody);
    end
end


