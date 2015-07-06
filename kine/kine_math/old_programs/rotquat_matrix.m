function [pf, q] = rotquat_matrix(model_3d, psi, theta, phi)
% Calculate quaternion components using the same angles that are input to
% rot3D
% psi = -psi;
%phi = -phi; Take this out to make frame rotate same as rot3D

e_0 = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2);
e_x = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2);
e_y = cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2);
e_z = cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2);

q = [e_0; e_x; e_y; e_z];
%-----------------------------------------
% Make rotation matrix from quaternion elements

Q = [ e_0^2+e_x^2-e_y^2-e_z^2 2*e_x*e_y-2*e_0*e_z     2*e_x*e_z+2*e_0*e_y;...
      2*e_x*e_y+2*e_0*e_z     e_0^2-e_x^2+e_y^2-e_z^2 2*e_y*e_z-2*e_0*e_x;...
      2*e_x*e_z-2*e_0*e_y     2*e_y*e_z+2*e_0*e_x     e_0^2-e_x^2-e_y^2+e_z^2 ];
  
pf = Q*model_3d;