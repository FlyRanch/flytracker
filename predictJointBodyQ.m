function [Pnew,tnew] = predictJointBodyQ(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function predicts the new body transformation using the twist
% representation
%
% the joint angles are predicted using the quaternion representation

persistent data dataB

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

% Joint indices
jidx = 7:12;

%global screw motion indices
gidx = 1:6;


P = Pfull(jidx,:);

%Convert joint angles to quaternions
Q = zeros(8,size(P,2));
for j = 1:size(P,2)
    LTheta = P(1:3,j);
    RTheta = P(4:6,j);
    Q(1:4,j) = quat2matNEW(Rx(LTheta(1))*Ry(LTheta(2))*Rz(LTheta(3)));
    Q(5:8,j) = quat2matNEW(Rx(RTheta(1))*Ry(RTheta(2))*Rz(RTheta(3)));
end
P = Q;

Pnew = zeros(size(Pfull,1),1);

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    load('exp083_ManualTrack.mat')
    dataB = [State(gidx,:) StartleState(gidx,:)];
    data = [StateQ(8:15,:) StartleStateQ(8:15,:)];
end


%this is the range of joint angle scalings
scale = 1;

tt = [1 size(data,2)];

minval = zeros(1,length(scale));
minidx = zeros(1,length(scale));

AllData = cell(length(scale),1);
AllTime = cell(length(scale),1);

for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %initialize
    scaleddata = zeros(size(data,1),length(scaledtime));
    %Scale the joint quaternion variables
    scaleddata(1:4,:) = slerpMAT(data(1:4,:),scale(k));
    scaleddata(5:8,:) = slerpMAT(data(5:8,:),scale(k));
    
    %scale the body twist data
    scaleddataB = zeros(size(dataB,1),length(scaledtime));
    for i = 1:size(dataB,1)
        scaleddataB(i,:) = interp1(tt(1):tt(2),dataB(i,:),scaledtime);
    end
    
    AllData{k} = [scaleddataB ;scaleddata];
    
    Error = zeros(1,length(scaledtime));
    
    for i = m:size(scaleddata,2)-1
        Error(i) = sum( sqrt( sum((scaleddata(:,i-m+1:i) - P).^2,1) ) );
    end

    [minval(k),minidx(k)] = min(Error(m:end-1));
    %fix the index because we ignore the first m-1 error values because
    %they are zero by construction.  Also ignore the last value in case it
    %is the prediction
    minidx(k) = minidx(k) + (m-1);
end

% scaleidx is the scale with the minimum error
[val,scaleidx] = min(minval);

% This is the n+1 state that is predicted
tsol = minidx(scaleidx);
%left wing
q1=AllData{scaleidx}(7:10,tsol);
q2=AllData{scaleidx}(7:10,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qL = quatprod(DeltaQ,P(1:4,end));


%right wing
q1=AllData{scaleidx}(11:14,tsol);
q2=AllData{scaleidx}(11:14,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qR = quatprod(DeltaQ,P(5:8,end));


JointPnew = [Rot2Joint(quat2matNEW(qL))
    Rot2Joint(quat2matNEW(qR))];

scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
    
%% Now calculate the predicted global transformation from prior data
% This is the global transform at the optimal time location
% Let's use Rosenhahn's notation

skew = @(p) [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];

P1twist = -AllData{scaleidx}(1:6,tsol);
P2twist = -AllData{scaleidx}(1:6,tsol+1);

%Relative motion in prior coordinate frame
T1 = screw2homo(P2twist)*inv(screw2homo(P1twist));
%T1 = inv(screw2homo(P1twist))*screw2homo(P2twist);
%Convert to twist
T1 = homo2screw(T1);

%Transformation from Prior coordinate frame to Current one
M1 = -Pfull(1:6,end);
G = screw2homo(M1)*inv(screw2homo(P1twist));
%G = inv(screw2homo(M1))*screw2homo(P1twist);

%Relative motion in current coordinate frame
T1prime = G * [skew(T1(4:6)) T1(1:3);zeros(1,4)] * inv(G);
%convert back to vector format
T1prime = [T1prime(1:3,4);T1prime(3,2);T1prime(1,3);T1prime(2,1)];
%scale the twist according to previous velocities
% the norm of the rotation axis defines the velocity
v_avg = mean( sqrt(sum(Pfull(4:6,:).^2,1)) );

%Calculate the full screw motion from origin in the current frame
%Use this to get the magnitude (i.e. velocity) for rescaling
M2 = inv(screw2homo(T1prime)*screw2homo(M1));
M2 = homo2screw(M2);

v_cur = norm(M2(4:6));

%keyboard
%rescale and transform again
T1prime = T1prime .* (v_avg/v_cur); 
%T1prime(4:6) = T1prime(4:6) .* (v_avg/v_cur); 

%M2 = screw2homo(T1prime)*screw2homo(M1);
%We store the forward kinematics
M2 = inv(screw2homo(T1prime)*screw2homo(M1));
%M2 = screw2homo(M1)*screw2homo(T1prime);
M2 = homo2screw(M2);
%keyboard
Pnew = [M2 ; JointPnew];


