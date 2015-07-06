%Plot other possible parameterizations for the fly wing joint.
%All of them have singularities.

%Euler angle rotation functions
Rx = @(theta) [1 0 0
    0 cos(theta) -sin(theta)
    0 sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 sin(theta)
    0 1 0
    -sin(theta) 0  cos(theta)];
Rz = @(theta) [cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];

load exp101_ManualTrack.mat


for k = 1:size(lw,2)
    lw = State(7:9,k);
    jj(:,k) = Rot2JointZXY(Rx(lw(1))*Ry(lw(2))*Rz(lw(3)));
end

figure; plot(jj(:,150:270)')
