function [Pnew,tnew] = predictMotion(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
global PAR
persistent data

P = Pfull(7:end,:);

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    %data = State(7:end,:);
    data = State;
end

% Joint indices
jidx = 7:12;
    
%this is the range of joint angle scalings
scale = 1.0;

tt = [1 size(data,2)];

minval = zeros(1,length(scale));
minidx = zeros(1,length(scale));

AllData = cell(length(scale),1);
for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %initialize
    scaleddata = zeros(size(data,1),length(scaledtime));
    for i = 1:size(data,1)
        scaleddata(i,:) = interp1(tt(1):tt(2),data(i,:),scaledtime);
    end
    
    AllData{k} = scaleddata;
    
    Error = zeros(1,length(scaledtime));
    
    %Don't include the last prior data values in case they are the the predicted ones 
    for i = m:size(scaleddata,2)-1
        Error(i) = sum( sqrt( sum((scaleddata(jidx,i-m+1:i) - P).^2,1) ) );
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
DeltaPnew = diff(AllData{scaleidx}(jidx,tsol:tsol+1),1,2);

%This is the predicted joint angles
if Error(minidx) > 0
    JointPnew = Pfull(jidx,end) + DeltaPnew;
else
    JointPnew = AllData{scaleidx}(jidx,tsol+1);
end
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
%keyboard


