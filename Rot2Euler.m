function Y = Rot2Euler(X)
%
% performs the inverse mapping from a rotation matrix or quaternion to the
% associated Euler angles, Rxyz.
% This inverse operation is taken from "Representing Attitude" J. Diebel 
% Tech. Report, pg. 12
%
% This Report is NON-STANDARD because it does not follow the right hand
% rule!!!! (see pg 5, Section 2.4)
%
% quaternion has form q = [q1 q2 q3 q0], where q0 is scalar part.
%
% Y = [phi theta psi]';

dim = numel(X);
switch dim
    case 9 %This is a Rotation Matrix
        Y(1,:) = atan2(X(2,3),X(3,3));
        Y(2,:) = -asin(X(1,3));
        Y(3,:) = atan2(X(1,2),X(1,1));
        
    case 4 %This is a quaternion
        q0 = X(4);
        q1 = X(1);
        q2 = X(2);
        q3 = X(3);
        
        Y(1,:) = atan2(2*(q2*q3+q0*q1),q3^2-q2^2-q1^2+q0^2);
        Y(2,:) = -asin(2*(q1*q3 - q0*q2));
        Y(3,:) = atan2(2*(q1*q2+q0*q3),q1^2+q0^2-q3^2-q2^2);
end