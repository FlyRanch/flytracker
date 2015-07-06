function [Pnew,tnew] = predictIndJoint(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function searches the prior database for each wing individually and
% not in a joint manner since the takeoff may have wing asymmetries that
% are better modeled separately
persistent data

% Joint indices
jidx = 7:12;
ljidx = 7:9;
rjidx = 10:12;
P = Pfull;

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    data = State;
end


%this is the range of joint angle scalings
scale = 0.5:0.1:2;

tt = [1 size(data,2)];

minvalLJ = zeros(1,length(scale));
minidxLJ = zeros(1,length(scale));

minvalRJ = zeros(1,length(scale));
minidxRJ = zeros(1,length(scale));

AllData = cell(length(scale),1);
for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %initialize
    scaleddata = zeros(size(data,1),length(scaledtime));
    for i = 1:size(data,1)
        scaleddata(i,:) = interp1(tt(1):tt(2),data(i,:),scaledtime);
    end
    
    AllData{k} = scaleddata;
    
    ErrorLJ = zeros(1,length(scaledtime));
    ErrorRJ = zeros(1,length(scaledtime));
    
    for i = m:size(scaleddata,2)-1
        ErrorLJ(i) = sum( sqrt( sum((scaleddata(ljidx,i-m+1:i) - P(ljidx,:)).^2,1) ) );
        ErrorRJ(i) = sum( sqrt( sum((scaleddata(rjidx,i-m+1:i) - P(rjidx,:)).^2,1) ) );
    end

    [minvalLJ(k),minidxLJ(k)] = min(ErrorLJ(m:end-1));
    [minvalRJ(k),minidxRJ(k)] = min(ErrorRJ(m:end-1));
    %fix the index because we ignore the first m-1 error values because
    %they are zero by construction.  Also ignore the last value in case it
    %is the prediction
    minidxLJ(k) = minidxLJ(k) + (m-1);
    minidxRJ(k) = minidxRJ(k) + (m-1);
end

tnew = [];
%% Left wing minimum
% scaleidx is the scale with the minimum error
[val,scaleidx] = min(minvalLJ);
% This is the n+1 state that is predicted
tsol = minidxLJ(scaleidx);
DeltaP_LJ = diff(AllData{scaleidx}(ljidx,tsol:tsol+1),1,2);

scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [tnew scaledtime(minidxLJ(scaleidx)) scale(scaleidx)];

%% Right wing minimum
% scaleidx is the scale with the minimum error
[val,scaleidx] = min(minvalRJ);
% This is the n+1 state that is predicted
tsol = minidxRJ(scaleidx);
DeltaP_RJ = diff(AllData{scaleidx}(rjidx,tsol:tsol+1),1,2);

Pnew = Pfull(:,end) + [zeros(6,1) ; DeltaP_LJ ; DeltaP_RJ];

scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [tnew scaledtime(minidxRJ(scaleidx)) scale(scaleidx)];
    



