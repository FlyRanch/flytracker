% Finds the norm of a quaternion
%
% q_in is a quaternion (4x1 or 1x4 vector)
% The output, N, is a scalar representing the norm (length) of the
% quaternion.

function N = normq(q_in)
N = sqrt(q_in(1)^2 + q_in(2)^2 + q_in(3)^2 + q_in(4)^2);
