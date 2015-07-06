% Quaternion multiplication
%
% A and B are both quaternions (1x4 or 4x1 vectors).
% The output, q, is also a quaternion.
% Note quaternion multiplication is not commutative, so multq(A,B) does not
% produce the same result as multq(B,A).

function q = multq(A,B)
o = 1; x = 2; y = 3; z = 4;
q(o) = A(o)*B(o) - A(x)*B(x) - A(y)*B(y) - A(z)*B(z);
q(x) = A(o)*B(x) + A(x)*B(o) + A(y)*B(z) - A(z)*B(y);
q(y) = A(o)*B(y) - A(x)*B(z) + A(y)*B(o) + A(z)*B(x);
q(z) = A(o)*B(z) + A(x)*B(y) - A(y)*B(x) + A(z)*B(o);
