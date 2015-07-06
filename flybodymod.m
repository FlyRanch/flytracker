function [x,y,z,s,th,X,Frenet,TTbody] = flybodymod(p,radius,len,PAR)
% [x,y,z,s,th,X,Frenet] = flybodymod(p,radius,len,PAR)
%
% 'p' - 2D B-spline points that define the bent centerline of the fly body
% in the X-Z plane.
%
% 'radius' - 1D B-spline points that define the body width along it's
% length.
%
% 'len' - the length of the centerline axis
%
% the function rotates the width profile around the centerline in the shape
% of a ellipse with axes [1 1.2];
%
% Note: I have rearranged the orientation of the generating curves to make
% it match with Gwyneth's body frame orientation

c = PAR.c; %4;
npts = size(p,1);
knotv = oknotP(length(p),c,len,PAR.etamax);

th = linspace(0,2*pi,PAR.T1);
% th = pi/2 is the ventral side of the fly (where B-vector is [0 0 1])
% th = 3*pi/2 is the dorsal side of the fly
s = linspace(-len/2,len/2,PAR.L1);

%Right now, assume that body centerline is not bent according to spline
%pts.

%Calculate Centerline and Frenet Frame
% [N,D1] = dbasisP(c,s,npts,knotv);
% X = N*p;

%Put the curve in the X-Z plane (head towards positive x, dorsal towards
%positive z)
% X = [X(:,2) zeros(size(X,1),1) X(:,1)];
% % 
% Tang = D1*p;
% Tang = [Tang(:,2) zeros(size(Tang,1),1) Tang(:,1)];
% MAG = sqrt(sum(Tang.^2,2));
% T = Tang./repmat(MAG,1,3);

X = [s' zeros(PAR.L1,1) zeros(PAR.L1,1)];
T = [ones(PAR.L1,1) zeros(PAR.L1,1) zeros(PAR.L1,1)];  

%Cross product with minus y axis
skew_my = -[0 0 1
    0 0 0
    -1 0 0];
%Cross product with z axis
skew_z = [0 -1 0
    1 0 0
    0 0 0];
N = (skew_my*T')';
B = repmat([0 -1 0],size(N,1),1);
%B = repmat([0 0 1],size(N,1),1);
R =  radiusspline1(s',radius,length(s),0,len);
%keyboard
TEMP = zeros(size(X));

for j = 1:length(th)
    for i = 1:3
        TEMP(:,i) = R.*(1.2*cos(th(j)).*B(:,i) + 1*sin(th(j)).*N(:,i));
    end
    RR = X + TEMP; 
    x(:,j) = RR(:,1);
    y(:,j) = RR(:,2);
    z(:,j) = RR(:,3);
end
%x = x + repmat(X(2,1),PAR.L1,length(th))
%drsl_idx ==> th = pi
lenidx = 2;
TTbody = X(lenidx,:) + R(lenidx)*(-cos(pi).*N(lenidx,:) + 1.2*sin(pi).*B(lenidx,:));

% X = X - repmat(TTbody,PAR.L1,1);
% x = x - repmat(TTbody(1),PAR.L1,length(th));
% y = y - repmat(TTbody(2),PAR.L1,length(th));
% z = z - repmat(TTbody(3),PAR.L1,length(th));

Frenet.T = T;
Frenet.N = N;
Frenet.B = B;
%keyboard