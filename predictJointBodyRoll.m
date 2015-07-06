function [Pnew,tnew] = predictJointBodyRoll(Pfull,N,PAR)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function predicts the new body transformation given the acclelration
% noise
%
% N is the acceleration noise given by the [a alpha]'

persistent data

% Joint indices
jidx = 7:12;

%global screw motion indices
gidx = 1:6;
%linear velocity indices
vidx = 13:15;
%angular velocity indices
omidx = 16:18;
%Roll velocity indices
ridx = 19;

%Noise value components
jointNoise = N(1:6);
a = N(7:9);
alpha = N(10:12);
dphi = N(13);

P = Pfull(jidx,:);
Pnew = zeros(size(Pfull,1),1);

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
Pold = Pfull(:,end);
%Add the new Joint locations
Pnew(jidx) = Pold(jidx) + DeltaPnew + jointNoise;

%Linear acceleration noise
Pnew(vidx) = Pold(vidx) + a*PAR.dt;

%angular acceleration noise
Pnew(omidx) = Pold(omidx) + alpha*PAR.dt;

%roll acceleration noise
Pnew(ridx) = Pold(ridx) + dphi*PAR.dt;

%amount of Roll
phi = Pold(ridx)*PAR.dt + 0.5*dphi*PAR.dt^2;
% Body frame 'x' axis
ttmp = screw2homo(Pold(gidx));
ttmp = ttmp(1:3,1);
%Calculate the tranformation 
G_roll = screw2homo([0 0 0 phi.*ttmp']);

%Calculate the linear displacement due to this acceleration noise
DT = Pold(vidx)*PAR.dt + 0.5*a*PAR.dt^2;

%Calculate the angular displacement due to  acceleration noise
DOmega = Pold(omidx)*PAR.dt + 0.5*alpha*PAR.dt^2;

%Calculate differential homogeneous transformation
%First compute rotation matrix using Rodriguez formula
theta = norm(DOmega);
if theta < eps
    % We will never have the pure translation case, so theta = 0 means
    % identity transformation
    R = eye(3);
else
    omega = DOmega/theta;
    
    beta = sin(theta);
    gamma = 1-cos(theta);
    omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        
    R = eye(3) + omegav*beta + omegav*omegav*gamma;
end

DG = [R DT;zeros(1,3) 1];

%apply transformations
Pnew(gidx) = homo2screw( DG*G_roll*screw2homo(Pold(gidx)) );

scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
    



