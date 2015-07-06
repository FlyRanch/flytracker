function [x,y,z] = flymod_BodyOnly(p,params,PAR)
% function [x,y,z] = flymod_BodyOnly(p,params,PAR)
%
% x,y,z are [1 1] cells that contain the body/head points, 
% 
% Twist notation is taken from Murray Book "Robotic Manipulation"
%
% p is the fly state denoted as
% p = [ xi          - global screw motion of the fly model xi = [v ; w]

%
% params is a structure that has all of the model parameters that encode
% the shape.

persistent xorg yorg zorg

x = cell(1,1);
y = x;
z = x;

p = p(:);
% the posture transformations
bodyscrew = p(1:6);
Gbody = screw2homo(bodyscrew);

Tbody = Gbody(1:3,4);
Rbody = Gbody(1:3,1:3);


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
    %headtilt = -18*(pi/180);
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


