% input: dataAngles = [psi(1)   psi(2)   psi(3)...
%                      theta(1) theta(2) theta(3)...
%                      phi(1)   phi(2)   phi(3)...]
%
% output: q = [e_0(1) e_0(2) e_0(3)...
%              e_x(1) e_x(2) e_x(3)...
%              e_y(1) e_y(2) e_y(3)...
%              e_z(1) e_z(2) e_z(3)...]

function quat = dataAngles2quat(dataAngles)

psi = dataAngles(1,:);
theta = dataAngles(2,:); 
phi = dataAngles(3,:);

e_0 = cos(phi/2).*cos(theta/2).*cos(psi/2) + sin(phi/2).*sin(theta/2).*sin(psi/2);
e_x = sin(phi/2).*cos(theta/2).*cos(psi/2) - cos(phi/2).*sin(theta/2).*sin(psi/2);
e_y = cos(phi/2).*sin(theta/2).*cos(psi/2) + sin(phi/2).*cos(theta/2).*sin(psi/2);
e_z = cos(phi/2).*cos(theta/2).*sin(psi/2) - sin(phi/2).*sin(theta/2).*cos(psi/2);


quat = [e_0; e_x; e_y; e_z];


