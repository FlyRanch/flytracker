function [Pnew,tnew] = predictMotionQNoise(Pfull,N,PAR)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function is for the quaternion representation of the state.
%
% It uses the scaled data to match and update the joint angles
%
% It uses the velocity estimates to integrate and update the body
% transformation
%
% This function predicts the new body transformation given the acclelration
% noise
% 
% See Julier 2007 "On Kalman filtering with nonlinear equality constraints"
% pg 2781 IEEE Trans on Sig. Process. for quaternion motion model.
%
% N is the acceleration noise given by the [a alpha]'

persistent data dataQ

% Joint indices
jidx = 8:15;

%global body transformation indices
gidx = 1:7;

P = Pfull(jidx,:);

Pnew = zeros(size(Pfull,1),1);

m = size(P,2);

if size(Pfull,1) > 15
    %linear velocity indices
    vidx = 16:18;
    %angular velocity indices
    omidx = 19:21;
else
    vidx = [];
    omidx = [];
end

%Noise value components
jointNoise = N(1:8);
a = N(9:11);
alpha = N(12:14);


% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    load('exp083_ManualTrack.mat')
    data = [StateQ(1:3,:) StartleStateQ(1:3,:)];
    dataQ = [StateQ(4:15,:) StartleStateQ(4:15,:)];
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
    scaleddataQ = zeros(size(dataQ,1),length(scaledtime));
    %Scale the body and joint quaternion variables
    scaleddataQ(1:4,:) = slerpMAT(dataQ(1:4,:),scale(k));
    scaleddataQ(5:8,:) = slerpMAT(dataQ(5:8,:),scale(k));
    scaleddataQ(9:12,:) = slerpMAT(dataQ(9:12,:),scale(k));
    
    %scale the body twist data
    scaleddata = zeros(size(data,1),length(scaledtime));
    for i = 1:size(data,1)
        scaleddata(i,:) = interp1(tt(1):tt(2),data(i,:),scaledtime);
    end
    
    AllData{k} = [scaleddata ;scaleddataQ];
    AllTime{k} = scaledtime;
    
    Error = zeros(1,length(scaledtime));
    
    for i = m:size(scaleddataQ,2)-1
        Error(i) = sum( sqrt( sum((scaleddataQ(5:12,i-m+1:i) - P).^2,1) ) );
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

%keyboard
%left wing orientation
q1=AllData{scaleidx}(8:11,tsol);
q2=AllData{scaleidx}(8:11,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qL = quatprod(DeltaQ,Pfull(8:11,end));


%right wing orientation
q1=AllData{scaleidx}(12:15,tsol);
q2=AllData{scaleidx}(12:15,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qR = quatprod(DeltaQ,Pfull(12:15,end));

Pold = Pfull(:,end);
%Add noise to the new Joint rotations
Pnew(jidx) = [qL ; qR] + jointNoise;

%Linear acceleration noise
Pnew(vidx) = Pold(vidx) + a*PAR.dt;

%angular acceleration noise
Pnew(omidx) = Pold(omidx) + alpha*PAR.dt;

%Calculate the linear displacement due to this acceleration noise
Pnew(1:3) = Pold(1:3) + Pold(vidx)*PAR.dt + 0.5*a*PAR.dt^2;

%Now, calculate the differential operator for the quaternion
aa = 0.5*PAR.dt*Pold(omidx(1)) + 0.25*PAR.dt^2*alpha(1);
bb = 0.5*PAR.dt*Pold(omidx(2)) + 0.25*PAR.dt^2*alpha(2);
cc = 0.5*PAR.dt*Pold(omidx(3)) + 0.25*PAR.dt^2*alpha(3);
dd = sqrt(aa^2+bb^2+cc^2);

D = eye(4).*cos(dd) + ...
    [0 cc -bb aa
    -cc 0 aa bb
    bb -aa 0 cc
    -aa -bb -cc 0].*(sin(dd)/dd);

Pnew(4:7) = D*Pold(4:7);
    
scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
