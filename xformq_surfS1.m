function [xx,yy,zz] = xformq_surfS1(x,y,z,T,q,scale)

%Perform coordinate transformation for a set of points that define a
%surface
if nargin == 4
    q = [1 0 0 0];
    scale = 1;
elseif nargin == 5
    scale = 1;
end

[nrow,ncol] = size(x);

%Make Row;
T = T(:)';

%Add zeros to work with quaternion math
pts = [zeros(nrow*ncol,1) reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];

% quaternion math uses formula from pg. 15 "Representing Attitude", by J. Diebel. 
% Quaternion is in the form [q0 q1 q2 q3] 
% where q0 is the scalar component.

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

Q = [q0 -q1 -q2 -q3
    q1 q0 q3 -q2
    q2 -q3 q0 q1
    q3 q2 -q1 q0];

Qconj = [q0 -q1 -q2 -q3
    q1 q0 -q3 q2
    q2 q3 q0 -q1
    q3 -q2 q1 q0];
%keyboard

% (Qconj'*Q * pts')' = pts*(Qconj'*Q)' matrix transpose
pts = pts*((Qconj'*Q)');
pts = scale.*pts(:,2:4) + repmat(T,nrow*ncol,1);

xx = reshape(pts(:,1),nrow,[]);
yy = reshape(pts(:,2),nrow,[]);
zz = reshape(pts(:,3),nrow,[]);