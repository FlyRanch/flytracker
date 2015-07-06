%convert joint state to quaternion state

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

for j = 1:size(State,2)
    p = State(:,j);

    bodyscrew = p(1:6);
%     Gbody = screw2homo(bodyscrew);
% 
%     Tbody = Gbody(1:3,4);
%     Rbody = Gbody(1:3,1:3);
%     q_body = quat2matNEW(Rbody);
% 
    LTheta = p(7:9);
    RTheta = p(10:12);

%     qL =quat2matNEW( Rx(LTheta(1))*Ry(LTheta(2))*Rz(LTheta(3)));
%     qR = quat2matNEW(Rx(RTheta(1))*Ry(RTheta(2))*Rz(RTheta(3)));
    qL = homo2screw([Rx(LTheta(1))*Ry(LTheta(2))*Rz(LTheta(3)) zeros(3,1);zeros(1,3) 1]);
    qR = homo2screw([Rx(RTheta(1))*Ry(RTheta(2))*Rz(RTheta(3)) zeros(3,1);zeros(1,3) 1]);


%     StartleStateQ(:,j) = [Tbody
%         q_body
%         qL
%         qR];
    nQ(:,j) = [bodyscrew
        qL(4:6)
        qR(4:6)];
%      nQ(:,j) = qL;

end