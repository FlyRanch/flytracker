function [Pnew] = predictQNoise(Pold,N,PAR)
%
%
% This function is for the quaternion representation of the state.
%
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

% Joint indices
jidx = 8:15;


Pnew = zeros(size(Pold,1),1);

if size(Pold,1) > 15
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


%Add noise to the new Joint rotations
Pnew(jidx) = Pold(jidx) + jointNoise;

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
if dd ~= 0
    D = eye(4).*cos(dd) + ...
        [0 cc -bb aa
        -cc 0 aa bb
        bb -aa 0 cc
        -aa -bb -cc 0].*(sin(dd)/dd);
else
    D = eye(4);
end

Pnew(4:7) = D*Pold(4:7);
    
