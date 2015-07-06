function [x,y,z] = flymod_BodyShapeOnly(p,X,params,PAR)
% function [x,y,z] = flymod(p,params,PAR)
%
% x,y,z are [1 3] cells that contain the body/head points, left wing, and
% right wing points, respectively;
% 
% Twist notation is taken from Murray Book "Robotic Manipulation"
%
% p is the fly state denoted as
% p = [ T_body      - global translation to center of fly body
%       q_body      - quaternion for body orientation
%       q_Lwing     - quaternion for left wing orientation
%       q_Rwing];   - quaternion for right wing orientation
%
%
% params is a structure that has all of the model parameters that encode
% the shape.

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

bodyrad = [0 ; X(1:18,:) ; 0];
headrad = [0 ; X(19:end,:) ; 0];
%=========================================================
%BODY
%=========================================================
%Instead of using the bent fly model from W. Dickson, let's just assume
%that the body axis is a straight line. Look inside 'flybodymod' for
%manipulation.
[xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(params.bodyctr,bodyrad,params.bodylen,PAR);

%=========================================================
%HEAD
%=========================================================
%Transformation to the neck;
% headtilt = 18*(pi/180);
headtilt = 0;
Ry = [cos(headtilt) 0 -sin(headtilt)
    0 1 0
    sin(headtilt) 0 cos(headtilt)];
Rneck = Ry*[Frenet.T(end,:)' Frenet.N(end,:)' Frenet.B(end,:)'];

%Tneck = X(end,:) + [-0.020 0 0]; %Shift head back a little
Tneck = X(end,:);

%Transformation to put reference frame at bottom of head and shift
% Thead = [params.headlen/2 0 -0.1];
Thead = [params.headlen/2 0 0];

[xhead,yhead,zhead] = flyheadmod(headrad,params.headlen,PAR);

%Translate reference frame at bottom of head
[xhead,yhead,zhead] = xform_surf(xhead,yhead,zhead,Thead);

%Transform to neck of body
[xhead,yhead,zhead] = xform_surf(xhead,yhead,zhead,Tneck,Rneck);

%Put head and body points together.
x{1} = [xbody ; xhead];
y{1} = [ybody ; yhead];
z{1} = [zbody ; zhead];


%=========================================================
%=========================================================

% global transformation for all body parts (and scaling!)
for k = 1:length(x)
    if k == 1 %body
        [x{k},y{k},z{k}] = xformq_surf(x{k},y{k},z{k},Tbody,qbody,params.bodyscale);
    else %wings
        [x{k},y{k},z{k}] = xformq_surf(x{k},y{k},z{k},Tbody,qbody,1);
    end
end


