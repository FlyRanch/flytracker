function [Pnew,tnew] = predictJointUpdateRoll(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
% 
% Only Predicts the wing joint variables.  
% The body transformation is updated based on the idea that the roll angle
% should be symmetric between the two wings (see 'fixRollTPlane.m').

global PAR 

persistent data

% Joint indices
jidx = 7:12;

P = Pfull(jidx,:);

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    data = State(jidx,:);
end


%this is the range of joint angle scalings
scale = 0.5:0.1:2.0;

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
DeltaPnew = diff(AllData{scaleidx}(:,tsol:tsol+1),1,2);

%keyboard
Pnew = Pfull(:,end) + [zeros(6,1) ; DeltaPnew];

%Roll Angle correction
Pnew = fixRollTPlane(Pnew,PAR);

scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
    



