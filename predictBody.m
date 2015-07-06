function Pnew = predictBody(Pfull,N,PAR)

% This function predicts the new body transformation given the acclelration
% noise
%
% N is the acceleration noise given by the [a alpha]'

% Joint indices
jidx = 7:12;

%global screw motion indices
gidx = 1:6;
%linear velocity indices
vidx = 13:15;
%angular velocity indices
omidx = 16:18;

%Noise value components
jointNoise = N(1:6);
a = N(7:9);
alpha = N(10:12);

Pnew = zeros(size(Pfull,1),1);


Pold = Pfull(:,end);
%Add the new Joint locations
Pnew(jidx) = Pold(jidx) + jointNoise;

%Linear acceleration noise
Pnew(vidx) = Pold(vidx) + a*PAR.dt;

%angular acceleration noise
Pnew(omidx) = Pold(omidx) + alpha*PAR.dt;

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

%apply transformation
Pnew(gidx) = homo2screw( DG*screw2homo(Pold(gidx)) );

    



