function Y = Rot2JointZXY(X)
%
% This assumes that the order of the wing rotations is reversed:
% Rz*Ry*Rx instead of Rx*Ry*Rz
%
% performs the inverse mapping from a rotation matrix or quaternion to the
% associated joint rotation Rxyz.
% This inverse operation is mathematica file "QuatTest.nb"
%
% quaternion has form q = [q1 q2 q3 q0], where q0 is scalar part.
%
% Y = [phi theta psi]';

dim = numel(X);
switch dim
    case 9 %This is a Rotation Matrix
        Y(1,:) = asin(X(3,2));
        Y(2,:) = atan2(-X(3,1),X(3,3));
        Y(3,:) = atan2(-X(1,2),X(2,2));
                
    case 4 %This is a quaternion
        error('Not working for quaternions yet!');
%         q0 = X(4);
%         q1 = X(1);
%         q2 = X(2);
%         q3 = X(3);
%         
%         Y(1,:) = atan2(2*(q2*q3+q0*q1),q3^2-q2^2-q1^2+q0^2);
%         Y(2,:) = -asin(2*(q1*q3 - q0*q2));
%         Y(3,:) = atan2(2*(q1*q2+q0*q3),q1^2+q0^2-q3^2-q2^2);
end