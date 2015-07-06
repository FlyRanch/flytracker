function Y = find_features1(x,N,lambda,C,thresh,IM,PHI)

%C - edge detection mask is row vector.
%lambda is a row vector of scalar multiples of the Normal vector

%These are the sample points along the normal vector direction for
%point x

%this uses a level set inside/outside image to determine if the
%edge point should correspond.

%Here we assume that lambda is >= 0 so we are only looking in the
%regions outside the worm

X = repmat(x,length(lambda),1) + repmat(lambda,2,1)'.*repmat(N,length(lambda), ...
                                                1);

%These are the corresponding image intensities for the sample points
Iout = getIntensity(X,IM);
%These are the corresponding intensity values for the body only segmented 
%images at the sample points.  
%Dout > 0 ==> inside
%Dout < 0 ==> outside
Dout = getIntensity(X,PHI);

%Convolve 1-D image with filter
F = imfilter(Iout,C','conv','replicate');
G = imfilter(Dout,C','conv','replicate');
idx = find(F > thresh);

[dum,bdymax] = max(G);
[dum,IMmax] = max(F);
if abs(IMmax-bdymax) < 5
    idx = [];
end
% for k = 1:length(idx)
%   %if the majority of points previous to the feature were outside
%   %then it is a valid match.
%   vec = 1:idx(k);
%   %ISvalid = round( length(find(Dout(vec) < 0)) / length(vec) );
%   
%   %if the last M points are not on the body, then it is a valid match
%   M = 3;  
% %   if ~isempty(Dout)
% %     ISvalid = (length(find(Dout(vec) > 0)) <= 2);
% %   else
% %     ISvalid = true;
% %   end
%   
%   %If there are any points previous to it that are inside the contour, then
%   %it is a valid match, otherwise, this worm is crossing over the other
%   %worm
%   if ~isempty(Dout)
%     ISvalid = (length(find(Dout(vec) > 0)) == 0);
%   else
%     ISvalid = true;
%   end
%   
%   if ISvalid
%     idx2keep = [idx2keep idx(k)];
%   end
% end

Y = X(idx,:);
