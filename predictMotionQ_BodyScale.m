function [Pnew,tnew] = predictMotionQ(Pfull)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function predicts the new body transformation using a translation
% vector and body quaternion
%
% the joint angles are predicted using the quaternion representation
%
% This functions make sure to ignore the jump
% discontinuities between the different prior motion patterns
persistent data dataQ transidx

% Joint indices
jidx = 8:15;

P = Pfull(jidx,:);

m = size(P,2);

% load data if it is not already in memory
if isempty(data)
    load('exp101_ManualTrack.mat')
    load('exp083_ManualTrack.mat')
    load('exp093_ManualTrack.mat')
    data = [StateQ_sym(1:3,:) StartleStateQ_sym(1:3,:) WideSweepStateQ(1:3,:)];
    dataQ = [StateQ_sym(4:15,:) StartleStateQ_sym(4:15,:) WideSweepStateQ(4:15,:)];
    transidx = [size(StateQ_sym,2) size(StartleStateQ_sym,2) size(WideSweepStateQ,2)];
    transidx = cumsum(transidx);
end


%this is the range of joint angle scalings
scale = [8/12 8/10 1];

tt = [1 size(data,2)];

minval = zeros(1,length(scale));
minidx = zeros(1,length(scale));

AllData = cell(length(scale),1);
AllTime = cell(length(scale),1);

for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %Determine the indices that are in the transition zone so we can ignore
    %them;
    spidx = [];
    for j = 1:length(transidx)
        spidx = [spidx find( transidx(j) <= scaledtime & ...
            scaledtime < transidx(j) + 1)];
    end
    
    
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
    
    % Artificial Fix - Make the Error in the transition zones arbitrarily
    % large so they can never get choosen.
    Error(spidx) = 10;
    
    %keyboard
    
    [minval(k),minidx(k)] = min(Error(m:end-1));
    %fix the index because we ignore the first m-1 error values because
    %they are zero by construction.  We also ignore the last one
    minidx(k) = minidx(k) + (m-1);
end

% scaleidx is the scale with the minimum error
[val,scaleidx] = min(minval);

% This is the n+1 state that is predicted
tsol = minidx(scaleidx);


%update body orientation
q1=AllData{scaleidx}(4:7,tsol);
q2=AllData{scaleidx}(4:7,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qB = quatprod(DeltaQ,Pfull(4:7,end));


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

%transformation
DT = diff(AllData{scaleidx}(1:3,tsol:tsol+1),1,2);
Tb = Pfull(1:3,end) + DT;


%% Scaled Body Motion
ttheta = 2*acos(qB(4));
omega = qB(1:3)./sin(ttheta/2);

beta = sin(ttheta);
gamma = 1-cos(ttheta);
omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];

R = eye(3) + omegav*beta + omegav*omegav*gamma;
A = (eye(3) - R)*omegav + omega*omega'.*ttheta;

v = inv(A)*Tb;

% Velocity from previous frames
theta = mean(2*acos(Pfull(7,:)));

qBnew = [omega.*sin(theta/2) ; cos(theta/2)];

beta = sin(theta);
gamma = 1-cos(theta);

R = eye(3) + omegav*beta + omegav*omegav*gamma;

%Now compute the Translation vector
Tnew = ( (eye(3) - R)*omegav + omega*omega'.*theta )*v;


Pnew = [Tnew
    qBnew
    qL
    qR];

keyboard
scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
    




