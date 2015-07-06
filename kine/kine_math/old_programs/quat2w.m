function [q_dot_seq, w_seq] = quat2wprime(q_seq)

% Find the derivative of the quaternion sequence
for j = 1:4
    sm_param = 0;
    y = q_seq(j,:);
    x = 1:size(y,2);

    [ sp, y_sm ] = spaps( x, y, sm_param );

    sp_der_1 = fnder( sp, 1 );
    y_der_1 = fnval(sp_der_1,x);

    q_dot_seq(j,:) = y_der_1;

end

% Calculate the angular velocity
for f = 1:size(exp.quat_dot,2)

    q = q_seq(:,f);
    q_dot = q_dot_seq(:,f);

    W = [-q(2)  q(1)  -q(4) q(3);...
        -q(3) q(4)  q(1)  -q(2);...
        -q(4)  -q(3) q(2)  q(1)];

    w_seq(:,f) = 2*W*q_dot;

end
