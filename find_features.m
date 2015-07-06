function Y = find_features(x,N,lambda,C,thresh,IM)

%C - edge detection mask is row vector.
%lambda is a row vector of scalar multiples of the Normal vector

%These are the sample points along the normal vector direction for
%point x
X = repmat(x,length(lambda),1) + repmat(lambda,2,1)'.*repmat(N,length(lambda),1);

%close all; imshow(IM); hold on
%plot(X(:,1),X(:,2),'ro')

%These are the corresponding image intensities for the sample points
Iout = getIntensity(X,IM);

%Convolve 1-D image with filter
F = imfilter(Iout,C','conv','replicate');
idx = find(F > thresh);

Y = X(idx,:);
%keyboard