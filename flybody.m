function [x,y,z,s,th,X,Frenet] = flybody(p,radius,len,PAR)
%[x,y,z,s,th,Y,Frenet] = modcurvespline(p,eta,wormradius,len)
%This approximates the bend angle of the worm with a bspline.  THe
%function parameters are the y coordinates of the bspline.

%warning off MATLAB:quadl:MinStepSize;
%warning off MATLAB:quadl:MaxFcnCount;


c = PAR.c; %4;
knotv = oknotP(length(p),c,len,PAR.etamax);

th = linspace(0,2*pi,10);
s = linspace(-len/2,len/2,PAR.L);

%Calculate Centerline and Frenet Frame
[N,D1,D2] = dbasisP(c,s,npts,knotv);
X = N*p;
X = [X zeros(size(X,1),1)];

Tang = D1*p;
Tang = [Tang zeros(size(Tang,1),1)];
MAG = sqrt(sum(Tang.^2,2));
T = Tang./repmat(MAG,1,3);

skew_z = [0 -1 0
          1 0 0  
          0 0 0];
      
N = (skew_z*T')';
B = repmat([0 0 1],size(N,1),1);

R =  radiusspline1(s',radius,length(s),0,len);
TEMP = zeros(size(X));

for j = 1:length(th)
    for i = 1:3
        TEMP(:,i) = R.*(-cos(th(j)).*N(:,i) + sin(th(j)).*B(:,i));
    end
    RR = X + TEMP;
    x(:,j) = RR(:,1);
    y(:,j) = RR(:,2);
    z(:,j) = RR(:,3);
end


Frenet.T = T;
Frenet.N = N;
Frenet.B = B;
