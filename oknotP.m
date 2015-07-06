function x = oknotP(n,k,len,etamax);

%Subroutine to generate a B-spline uniform periodic knot vector 
%    k            = order of the basis function (degree + 1,
%    i.e. cubic => k=4)
%    n            = the number of defining polygon vertices
%    len          = length of the worm

%    x            = knot vector


x = zeros(1,n+k);

% I want the knot vector to be such that is satisfies the
% constraint of xmin <= t < xmax.       (1)

% curve parameter is t = linspace(0,LengthOfWorm), so if we want to
% be able to shift the curve along it's parameter length, and
% calculate P(t + eta) using the same control points, then the knot
% vector, x, must be so that t+eta always satisfies (1) whether x
% is < or > 0. 
%So, we want it to go from [-etamax len+etamax];
% OR
% we want it to go from [-etamax-len/2 len/2+etamax];
epsilon = 0;

% the total number of points in the knot vector is n+k
% for a periodic b-spline, with knot vector xx = [0 1 2 3 4 5 6 7] 
% (n = 4, k = 4), the parameter range is 3 <= t <= 4, i.e the first
% and last k-1 knot spans are not used.

% I want the parameter range to go from -len/2:len/2 and it should
% contain (n+k) - 2*(k-1) values uniformly spaced.
%if etamax >=0
  mid = linspace(-etamax-len/2,etamax+len/2,n-k+2);

%else
%  mid = linspace(etamax-len/2,len/2,n-k+2);
%end

delta = len / (n-k+2 - 1);
beg = mid(1) + (-1*(k-1):-1)*delta ;
endd = mid(end) + (1:(k-1))*delta;


x = [beg mid endd];
