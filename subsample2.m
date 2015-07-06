function A = subsample2(I,n)

% A = subsample(I,n) subsamples an image B by a factor of '2^n' where
% 'n' is an integer.
% 
if ~isa(I,'double')
    I = double(I);
end

Arow = [];
Acol = [];

%check that subsample is evenly divisible
%if (mod(size(I,1),n) ~= 0 & n ~=0)
%  error('Subsample factor is not evenly divisible into rows');
%end

A = I;


t = 1;
tt = 2;

for k = 1:n
  A = imfilter(A,fspecial('gaussian',[3 3],1),'replicate','conv');
  B = A(t:tt:end,t:tt:end);
  A = B;
end

