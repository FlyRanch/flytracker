function [x,y,z,s,th,X,Frenet] = flyheadmod(radius,len,PAR)
% [x,y,z,s,th,X,Frenet] = flyheadmod(radius,len,PAR)
%
% Takes the outline profile defined by the spline points in 'radius' 
% (in the X-Z plane) and rotates them around the X-axis according to an 
% ellipse with axes [1 1.2] to create a mushroom shaped head.
% 
% 'radius' 2D spline control points
% 'len' is the height of the head.
%
% Note: I have rearranged the orientation of the generating curves to make
% it match with Gwyneth's body frame orientation
th = linspace(0,2*pi,PAR.T1);
s = linspace(-len/2,len/2,PAR.L2);

%Calculate Centerline and Frenet Frame
X = [s' zeros(PAR.L2,1) zeros(PAR.L2,1)];
T = [ones(PAR.L2,1) zeros(PAR.L2,1) zeros(PAR.L2,1)];      

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
TEMP = zeros(size(T));

for j = 1:length(th)
    for i = 1:3
        TEMP(:,i) = R.*(1.2*cos(th(j)).*B(:,i) + 1*sin(th(j)).*N(:,i));
    end
    RR = X + TEMP;
    x(:,j) = RR(:,1);
    y(:,j) = RR(:,2);
    z(:,j) = RR(:,3);
end


Frenet.T = T;
Frenet.N = N;
Frenet.B = B;
