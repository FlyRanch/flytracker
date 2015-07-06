function x = oknot(n,k,len,etamax);

%Subroutine to generate a B-spline uniform open knot vector with
%multiplicity equal to the order at the ends.
%    k            = order of the basis function (degree + 1,
%    i.e. cubic => k=4)
%    n            = the number of defining polygon vertices
%    x            = knot vector


x = zeros(1,n+k);
% $$$ x(1) = 0;
% $$$ for i = 2:n+k
% $$$   if ( (i > k) & (i < n+2))
% $$$     x(i) = x(i-1) + 1;
% $$$   else
% $$$     x(i) = x(i-1);
% $$$   end
% $$$ end
% $$$ 
% $$$ %normalize the knot vector
% $$$ x = x./max(x);

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
% $$$ if etamax < 0
% $$$   beg = repmat(etamax,1,k-1);
% $$$   endd = repmat(len,1,k-1);
% $$$ elseif etamax > 0
% $$$   beg = repmat(0,1,k-1);
% $$$   endd = repmat(len+etamax,1,k-1);
% $$$ else
% $$$   beg = repmat(0,1,k-1);
% $$$   endd = repmat(len,1,k-1);
% $$$ end

beg = repmat(-etamax-len/2,1,k-1);
endd = repmat(len/2+etamax,1,k-1);

x = [beg linspace(beg(1),endd(1),n-k+2) endd];
