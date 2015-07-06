% Rotate a body in 3D space using a quaternion
%
% q is a unit quaternion (4x1 or 1x4)
% origVector is the object to be rotated (3xn)

function rotVector=rot3D(origVector,q)

Q = q2Q(q);
rotVector = Q*origVector;

%% OLD Kine Rotation Scheme
% ps = theta;
% th = phi;
% ph = alpha;
% 
% %Rotate the wing around x axis.
% m1=[[1 0 0];[0 cos(ph) -sin(ph)];[0 sin(ph) cos(ph)]];
% v2 = m1 * origVector;
% 
% %Rotate the wing around y axis.
% m2=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
% v3 = m2 * v2;
% 
% %Rotate the wing around z axis.
% m3=[[cos(ps) -sin(ps) 0];[sin(ps) cos(ps) 0];[0 0 1]];
% v4 = m3 * v3;
% 
% transVector=v4;

%% OLD OLD Kine Rotation Scheme

% First rotation (theta) is around z clockwise + (yielding x', y', z'=z)
% Second rotation (phi) is around y' ANTIclockwise + (yielding x'', y''=y', z'')
% Third rotation (alpha) is around x'' ANTIclockwise + (yielding x'''=x'', y''', z''')

% %Rotate the wing around x axis.
% m1=[[1 0 0];[0 cos(alpha) -sin(alpha)];[0 sin(alpha) cos(alpha)]];
% v2 = m1 * origVector;
% 
% %Rotate the wing around y axis.
% m2=[[cos(phi) 0 sin(phi)];[0 1 0];[-sin(phi) 0 cos(phi)]];
% v3 = m2 * v2;
% 
% %Rotate the wing around z axis.
% m3=[[cos(theta) sin(theta) 0];[-sin(theta) cos(theta) 0];[0 0 1]];
% v4 = m3 * v3;

%transVector=v4;
