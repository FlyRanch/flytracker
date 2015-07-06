function G = screw2homo(xi)

% Convert a screw motion to the associated homogeneous transformation
%
% xi = [v ; w]
%
% G = [R T
%      0 1];
%
%

xi = xi(:); %make column
v = xi(1:3);
w = xi(4:6);

bigeps = 10e+20*eps;

%First compute rotation matrix using Rodriguez formula
theta = norm(w);
if theta < eps
    % We will never have the pure translation case, so theta = 0 means
    % identity transformation
    R = eye(3);
    T = zeros(3,1);
else
    omega = w/theta;
    
    beta = sin(theta);
    gamma = 1-cos(theta);
    omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        
    R = eye(3) + omegav*beta + omegav*omegav*gamma;
    
    %Now compute the Translation vector
    T = ( (eye(3) - R)*omegav + omega*omega'.*theta )*v;

end

G = [R T
    zeros(1,3) 1];
