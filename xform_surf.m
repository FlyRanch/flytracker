function [xx,yy,zz] = xform_surf(x,y,z,T,R,scale)

%Perform coordinate transformation for a set of points that define a
%surface
%
% Could be written better with homogeneous representation

if nargin == 4
    R = eye(3);
    scale = 1;
elseif nargin == 5
    scale = 1;
end

[nrow,ncol] = size(x);

%Make Row;
T = T(:)';

pts = [reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];

pts = scale.*(pts*R') + repmat(T,nrow*ncol,1);
xx = reshape(pts(:,1),nrow,[]);
yy = reshape(pts(:,2),nrow,[]);
zz = reshape(pts(:,3),nrow,[]);