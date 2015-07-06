function Q = quatprod(q1,q2)

%This function will perform quaternion multiplication for 2
%quaternions q1 * q2, .  Remember, multiplication is not commutative
% Q = quatprod(q1,q2)

%Make v2 col
q2 = q2(:);

q1plus = [q1(4) -q1(3) q1(2) q1(1);
    q1(3) q1(4) -q1(1) q1(2);
    -q1(2) q1(1) q1(4) q1(3);
    -q1(1) -q1(2) -q1(3) q1(4)];

%make column
Q = q1plus*q2;


