function Iout = getIntensity(X,im)

%Iout = getIntensity(X,im)
%This function calculates the intensity values in 'im' and the
%points X;  X in a [N 2] matrix of points and it assumes that the
%coordinate axis is in the lower left. 

%Iout is a [N 1] vector of intensity values at the corresponding
%points based on bilinear interpolation.

[nrow,ncol] = size(im);

%Convert X into (row,col) coordinates
Xidx = [X(:,2) X(:,1)];

%Find the indices of the four neighboring pixels by Rounding
%[rmin rmax cmin cmax]

IDX = [floor(Xidx(:,1)) ceil(Xidx(:,1)) floor(Xidx(:,2)) ceil(Xidx(:,2))];

%Get rid of sample points that are outside the image boundary 
ii = find(IDX(:,1) <= 0);
IDX(ii,:) = [];
Xidx(ii,:) = [];

ii = find(IDX(:,2) > nrow);
IDX(ii,:) = [];
Xidx(ii,:) = [];

ii = find(IDX(:,3) <= 0);
IDX(ii,:) = [];
Xidx(ii,:) = [];

ii = find(IDX(:,4) > ncol);
IDX(ii,:) = [];
Xidx(ii,:) = [];

%Calculate matrix of weights w = (1-abs(x-i))*(1-abs(y-j))
% W is a [N 4] matrix
W = [(1-abs(Xidx(:,1) - IDX(:,1))).*(1-abs(Xidx(:,2) - IDX(:,3)))...
     (1-abs(Xidx(:,1) - IDX(:,1))).*(1-abs(Xidx(:,2) - IDX(:,4)))...
     (1-abs(Xidx(:,1) - IDX(:,2))).*(1-abs(Xidx(:,2) - IDX(:,3)))...
     (1-abs(Xidx(:,1) - IDX(:,2))).*(1-abs(Xidx(:,2) - IDX(:,4)))];

%now create a matrix of image intensities with dimension [N 4] to
%perform the inner product

I = zeros(size(W));
for k = 1:size(I,1)
  I(k,:) = [im(IDX(k,1),IDX(k,3)) im(IDX(k,1),IDX(k,4)) im(IDX(k,2),IDX(k,3)) ...
            im(IDX(k,2),IDX(k,4))];
end

Iout = sum(W.*I,2);
