% Rotate a body in 3D space. This function has the opposite effect of rot3D while taking the same input arguments. 
% Applying rot3D and subsequently invRot3D with the same angle arguments produces the original vector. 

function transVector=invRot3D(origVector,theta,phi,alpha)

%alpha = pi*1/8;
%phi   =-pi*1/4;
%theta = pi*4/8;

%Rotate the wing around z axis.
m1=[[cos(theta) sin(theta) 0];[-sin(theta) cos(theta) 0];[0 0 1]];
v2 = m1' * origVector;

%Rotate the wing around y axis.
m2=[[cos(phi) 0 sin(phi)];[0 1 0];[-sin(phi) 0 cos(phi)]];
v3 = m2' * v2;

%Rotate the wing around x axis.
m3=[[1 0 0];[0 cos(alpha) -sin(alpha)];[0 sin(alpha) cos(alpha)]];
v4 = m3' * v3;

transVector=v4;


