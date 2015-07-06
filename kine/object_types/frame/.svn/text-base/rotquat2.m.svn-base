function [pf] = rotquat2(p, psi, theta, phi)

% NOTE: THIS SCHEME IS INCORRECT AND NO LONGER WORKS WHEN USE THE CORRECT
% INPUTS!!!!!!

% Calculate quaternion components using the same angles that are input to
% rot3D
%psi = -psi;
%phi = -phi; Take this out to make frame rotate same as rot3D

e_0 = cos(phi/2)*cos(theta/2)*cos(psi/2) - sin(phi/2)*sin(theta/2).*sin(psi/2);
e_x = sin(phi/2)*cos(theta/2)*cos(psi/2) + cos(phi/2)*sin(theta/2).*sin(psi/2);
e_y = cos(phi/2)*sin(theta/2)*cos(psi/2) - sin(phi/2)*cos(theta/2).*sin(psi/2);
e_z = -cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2).*cos(psi/2);

%-----------------------------------------
% Rotate using quaternion
%p' = (q)(p)(q_inv)
q = [e_0 e_x e_y e_z]';

pf = multq(multq(q,p),conjq(q))';

%---------------------------------------------------------
function N = normq(q_in)
N = sqrt(q_in(1)^2 + q_in(2)^2 + q_in(3)^2 + q_in(4)^2);

function q = multq(A,B)
o = 1; x = 2; y = 3; z = 4;
q(o) = A(o)*B(o) - A(x)*B(x) - A(y)*B(y) - A(z)*B(z);
q(x) = A(o)*B(x) + A(x)*B(o) + A(y)*B(z) - A(z)*B(y);
q(y) = A(o)*B(y) - A(x)*B(z) + A(y)*B(o) + A(z)*B(x);
q(z) = A(o)*B(z) + A(x)*B(y) - A(y)*B(x) + A(z)*B(o);

function q_out = conjq(q_in)
q_out = [q_in(1) -q_in(2) -q_in(3) -q_in(4)];