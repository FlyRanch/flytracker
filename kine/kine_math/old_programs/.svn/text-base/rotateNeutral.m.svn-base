% Function rotateNeutral
% 
% CALLING FUNCTION: 
% ACTIONS: Takes vector in axis and "un"-rotates it so the head is at (0,1)
% and the tail at (0,0)
% INPUTS: body_axis = [ x1 y1 z1; x2 y2 z2 ]
%                   = [ head_coords; tail_coords]
% PARENT PROGRAM: Kine_v2_1
% LAST MODIFIED: February 10, 2004 by gwyneth


function data_pts_rot = rotateNeutral(body_axis)

axis_v = body_axis(2,:) - body_axis(1,:);
[th,phi,r] = cart2sph(axis_v(1), axis_v(2), axis_v(3));




% Rotate a body in 3D space. 

% First rotation (theta) is around z clockwise + (yielding x', y', z'=z)
% Second rotation (phi) is around y' ANTIclockwise + (yielding x'', y''=y', z'')
% Third rotation (alpha) is around x'' ANTIclockwise + (yielding x'''=x'', y''', z''')


function transVector=rot3D(origVector,theta,phi,alpha)

%alpha = pi*1/8;
%phi   =-pi*1/4;
%theta = pi*4/8;


%Rotate the wing around x axis.
m1=[[1 0 0];[0 cos(alpha) -sin(alpha)];[0 sin(alpha) cos(alpha)]];
v2 = m1 * origVector;

%Rotate the wing around y axis.
m2=[[cos(phi) 0 sin(phi)];[0 1 0];[-sin(phi) 0 cos(phi)]];
v3 = m2 * v2;

%Rotate the wing around z axis.
m3=[[cos(theta) sin(theta) 0];[-sin(theta) cos(theta) 0];[0 0 1]];
v4 = m3 * v3;


transVector=v4;
