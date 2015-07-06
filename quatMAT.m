function Q = quatMAT(q1,tag)

%This function will perform quaternion multiplication for 2
%quaternions q1 * q2, .  Remember, multiplication is not commutative
% Q = quatprod(q1,q2)

if nargin < 2
    
    Q = [q1(4) -q1(3) q1(2) q1(1);
        q1(3) q1(4) -q1(1) q1(2);
        -q1(2) q1(1) q1(4) q1(3);
        -q1(1) -q1(2) -q1(3) q1(4)];
    
elseif strcmp(tag,'r')
    
    Q = [q1(4) q1(3) -q1(2) q1(1);
       -q1(3) q1(4) q1(1) q1(2);
        q1(2) -q1(1) q1(4) q1(3);
        -q1(1) -q1(2) -q1(3) q1(4)];
end
    





