function [x,y,z] = flymod1(p,ttemp,params,PAR)
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


% PAR.dt = 1/6000;
% PAR.etamax = 0;
% PAR.c = 4;
% PAR.pixpermm = 1;
% PAR.numfly = 1;
% PAR.L = 40;

x = cell(3,1);
y = x;
z = x;

%=========================================================
%BODY
%=========================================================
[xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(params.bodyctr,params.bodyrad,params.bodylen,PAR);

%=========================================================
%HEAD
%=========================================================
%Transformation to the neck;
headtilt = 18*(pi/180);
Ry = [cos(headtilt) 0 -sin(headtilt)
    0 1 0
    sin(headtilt) 0 cos(headtilt)];
Rneck = Ry*[Frenet.T(end,:)' -Frenet.B(end,:)' Frenet.N(end,:)'];

Tneck = X(end,:) + [-0.020 0 0]; %Shift head back a little

%Transformation to put reference frame at bottom of head and shift
Thead = [params.headlen/2 0 -0.1];

[xhead,yhead,zhead] = flyheadmod(params.headrad,params.headlen,PAR);

%Translate reference frame at bottom of head
[xhead,yhead,zhead] = xform_surf(xhead,yhead,zhead,Thead);

%Transform to neck of body
[xhead,yhead,zhead] = xform_surf(xhead,yhead,zhead,Tneck,Rneck);

%=========================================================
%WING
%=========================================================
[x{2},y{2},z{2}] = flywingmod(params.wing,params.winglen,PAR);

%Make Right Wing too
x{3} = x{2};
y{3} = y{2};
z{3} = z{2};

%Approximate translation from center of fly body to the wing joint location
%JTrans = [0 +/-0.35 0.5]
%RJTrans = [-0.18 -0.37 0.52];
%LJTrans = [-0.18 0.37 0.52];
RJTrans = [0 0 0];
LJTrans = [0 0 0];

%Rotate about z-axis by 180 degrees
Rz = @(theta) [cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];

% Translate Left wing
[x{2},y{2},z{2}] = xform_surf(x{2},y{2},z{2},LJTrans);

% Translate Right wing
[x{3},y{3},z{3}] = xform_surf(x{3},y{3},z{3},RJTrans);

%=========================================================
%=========================================================
%Put head and body points together.
x{1} = [xbody ; xhead];
y{1} = [ybody ; yhead];
z{1} = [zbody ; zhead];

% Now apply the posture transformations
Tbody = p(1:3);
qbody = p(4:7);
qLwing = p(8:11);
qRwing = p(12:15);

T_L = ttemp(1:3);
T_R = ttemp(4:6);

% global transformation for all body parts (and scaling!)
for k = 1:length(x)
    if k == 1
        [x{k},y{k},z{k}] = xformq_surf(x{k},y{k},z{k},Tbody,qbody,params.bodyscale);
    else
        %[x{k},y{k},z{k}] = xformq_surf(x{k},y{k},z{k},Tbody,qbody);
    end
end

% Now apply wing transformations
[x{2},y{2},z{2}] = xformq_surf(x{2},y{2},z{2},T_L,qLwing,params.wingscale);
[x{3},y{3},z{3}] = xformq_surf(x{3},y{3},z{3},T_R,qRwing,params.wingscale);
