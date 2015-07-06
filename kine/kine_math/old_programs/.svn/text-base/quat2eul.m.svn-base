% Convert quaternion to Euler angles
%
% Input, q, is a unit quaternion representing a rotation.
% The outputs are Euler angles in radians representing a rotation around
% axes in the lab frame, in the following order:
%   1) bank (rotation around x-axis)
%   2) elev (rotation around y-axis)
%   3) head (rotation around z-axis)
%
% This corresponds to a rotation matrix composed as: RzRyRx
%
% Note: Euler angles can be thought of as representing rotations around
% fixed axes in the lab frame OR as representing rotations around
% non-orthogonal axes that rotate sequentially as the object is rotated.

function [bank, elev, head] = quat2eul(q)

e0 = q(1);
ex = q(2);
ey = q(3);
ez = q(4);

% Quaternion to Euler
if round( 10^9*(e0*ey-ex*ez))/10^9 == 0.5
    
    elev = pi/2;
    head = 0;
    bank = 2*asin(ex/cos(pi/4)) + head;
    
elseif round( 10^9*(e0*ey-ex*ez))/10^9 == -0.5
    
    elev = -pi/2;
    head = 0;
    bank = 2*asin(ex/cos(pi/4)) - head;
    
else
    
    elev = asin( 2*(e0*ey-ex*ez) );
    head = atan2( 2*(e0*ez + ex*ey), e0^2+ex^2-ey^2-ez^2 );
    bank = atan2( 2*(e0*ex + ey*ez), e0^2+ez^2-ex^2-ey^2 );
    
end