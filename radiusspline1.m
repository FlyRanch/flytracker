function R = radiusspline1(ss,B,N,eta,len)

% 
% d is the maximum radius from centerline

% Create a vector of radius values based on size of s, and pick off the
% value that is closest to where we are along the centerline.  Radius
% values taper to zero at the ends.

npts = length(B);
c = 4;
knotv = oknot(npts,c,len,eta);
[N,D1,D2] = dbasis(c,ss,npts,knotv);

R = zeros(size(ss));

%Interpolate the radius value;
R = N*B;

