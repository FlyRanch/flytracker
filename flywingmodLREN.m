function [x,y,z,s,v,X,Frenet] = flywingmodLREN(p,len,PAR)
% [x,y,z,s,th,X,Frenet] = flywingmodLREN(p,len,PAR)
%
% 'p' - B-spline control points that define a closed curve that is the
% profile of wing
%
% 'len' is the thickness of the wing, but I currently assume it is planar
%
% the function takes the profile and scales it down to zero (the center of
% the wing
%
% Note: I have rearranged the orientation of the generating curves to make
% it match with Gwyneth's wing model orientation

c = PAR.c; %4;
npts = size(p,1);
knotv = oknotP(length(p),c,1,PAR.etamax);

%Make the wing 4 times thicker to have is show in pixels.
thickness = 0.0015;
w = linspace(-1,1,3);
v = linspace(1,0,PAR.T2);
s = linspace(-.5,.5,PAR.L3);

%Calculate wing boundary curve
[N,D1,D2] = dbasisP(c,s,npts,knotv);
X = N*p;

%Put the curve in the X-Y plane
X = [X zeros(size(X,1),1)];

%=================================================
%Change the orientation to match the left side of the fly
%x => -y
%y => x
%z => -z
X = [X(:,2) -X(:,1) X(:,3)];
%=================================================

Xp = X;


Tang = D1*p;
Tang = [Tang(:,1) Tang(:,2) zeros(size(Tang,1),1)];
MAG = sqrt(sum(Tang.^2,2));
T = Tang./repmat(MAG,1,3);

skew_z = [0 -1 0
          1 0 0  
          0 0 0];
      
N = (skew_z*T')';
Np = N(:,2:3);
B = repmat([0 0 1],size(N,1),1);

%Calculate curvature from dot product of second derivative and
%normal vector MAG(s) * Tang(s) = D1*B
%              MAG(s) * Kappa(s) * Nrml(s) = D2*B

% kappa = sum( (D2*p) .* Np,2) ./MAG;
% kk = [kappa kappa];
% TEMP = zeros(size(X));

% shift all points so that the origin is located at the joint 
% The joint of the fly wing is located at [1.294 -0.1225 0]
%wingTrans = repmat([1.294 -0.1225 0],PAR.L3,1);
wingTrans = repmat([-0.1225 -1.294 0],PAR.L3,1);
Zshift = repmat([0 0 thickness],PAR.L3,1);
for j = 1:length(w)
    RR = Xp + wingTrans + w(j).*Zshift;
    x(:,j) = RR(:,1);
    y(:,j) = RR(:,2);
    z(:,j) = RR(:,3);
end


Frenet.T = T;
Frenet.N = N;
Frenet.B = B;
