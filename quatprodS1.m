function qp = quatprodS1(q,p)

%This function will perform quaternion multiplication for 2
%quaternions q1 * q2, .  Remember, multiplication is not commutative
% Q = quatprod(q1,q2)
%
% Rearrange so that scalar part of quaternion is 1st element
% See page 14 of Diebel Technical Report
%Make v2 col
p = p(:);

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

Q = [q0 -q1 -q2 -q3
    q1 q0 q3 -q2
    q2 -q3 q0 q1
    q3 q2 -q1 q0];

%make column
qp = Q*p;


