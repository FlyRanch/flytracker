function [Pnew,tnew] = predictJointLoc(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function also predicts the wing joint locations based on the prior
% database
persistent data

% Joint angle indices
jidx = 7:12;

% Joint location indices
jlidx = 13:18;

P = Pfull(jidx,:);

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    data = [State ; JointDisplacement];
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

% This is the n+1 state that is predicted (joint angle and joint location
% displacement)
tsol = minidx(scaleidx);
DeltaPnew = diff(AllData{scaleidx}(jidx,tsol:tsol+1),1,2);

%keyboard
Pnew = Pfull(:,end) + [zeros(6,1) ; DeltaPnew ; zeros(6,1)];

%Replace the joint locations with the average of the predicted one from 
%the database and the current estimate.
Pnew(jlidx) = AllData{scaleidx}(jlidx,tsol+1);
%Pnew(jlidx) = mean([Pnew(jlidx) AllData{scaleidx}(jlidx,tsol+1)],2);

scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
    



