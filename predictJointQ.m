function [Pnew,tnew] = predictJointQ(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
% 
% Only Predicts the wing joint variables.  The body transformation is left
% the same
%
% This function transforms the joint variables into a quaternion, performs
% the prediction in quaternion space (no discontinuities), and then
% transforms back to joint angles


persistent data

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

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    load('exp083_ManualTrack.mat')
    data = [StateQ(8:15,:) StartleStateQ(8:15,:)];
end


%this is the range of joint angle scalings
scale = 1;
%scale = 1;
tt = [1 size(data,2)];

minval = zeros(1,length(scale));
minidx = zeros(1,length(scale));


AllData = cell(length(scale),1);
AllTime = cell(length(scale),1);
for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %initialize
    scaleddata = zeros(size(data,1),length(scaledtime));
    scaleddata(1:4,:) = slerpMAT(data(1:4,:),scale(k));
    scaleddata(5:8,:) = slerpMAT(data(5:8,:),scale(k));
    
    
    AllData{k} = scaleddata;
    AllTime{k} = scaledtime;
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
q1=AllData{scaleidx}(1:4,tsol);
q2=AllData{scaleidx}(1:4,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qL = quatprod(DeltaQ,P(1:4,end));


%right wing
q1=AllData{scaleidx}(5:8,tsol);
q2=AllData{scaleidx}(5:8,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qR = quatprod(DeltaQ,P(5:8,end));


Pnew = [Pfull(1:6,end) 
    Rot2Joint(quat2matNEW(qL))
    Rot2Joint(quat2matNEW(qR))];

scaledtime = AllTime{scaleidx};
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
%keyboard




