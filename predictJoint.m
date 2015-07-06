function [Pnew,tnew] = predictJoint(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
% 
% Only Predicts the wing joint variables.  The body transformation is left
% the same

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
%scale = 1;
tt = [1 size(data,2)];

minval = zeros(1,length(scale));
minidx = zeros(1,length(scale));

%Determine the time indices where the joint angles have a discontunity (at
%+- pi from atan2)
didx = cell(6,1);
Dtheta = diff(data,1,2);
for i = 1:size(Dtheta,1)
    didx{i} = find(abs(Dtheta(i,:)) >= pi);
end


AllData = cell(length(scale),1);
AllTime = cell(length(scale),1);
for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %initialize
    scaleddata = zeros(size(data,1),length(scaledtime));
    for i = 1:size(data,1)
        scaleddata(i,:) = interp1(tt(1):tt(2),data(i,:),scaledtime);
    end
    
    
    % Now, get rid of the indices that are in the region of discontinuity
    idx2remove = [];
    for i = 1:size(scaleddata,1)
        if ~isempty(didx{i})
            for j = 1:length(didx{i})
                beg = didx{i}(j);
                tmpidx = find(beg < scaledtime & scaledtime < beg+1);
                %Set these values that are in the region of discontinuity
                %to the value that is just outside the discontinuity.
                if ~isempty(tmpidx)
                    %scaleddata(i,tmpidx) = scaleddata(i,tmpidx(end)+1);
                    scaleddata(i,tmpidx) = data(i,beg+1);
                end
            end
        end
    end
    
    %scaleddata(:,idx2remove) = [];
    %scaledtime(:,idx2remove) = [];
    
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
DeltaPnew = diff(AllData{scaleidx}(:,tsol:tsol+1),1,2);


Pnew = Pfull(:,end) + [zeros(6,1) ; DeltaPnew];

scaledtime = AllTime{scaleidx};
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];

keyboard



