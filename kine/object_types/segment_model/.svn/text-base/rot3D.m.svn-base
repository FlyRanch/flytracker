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
