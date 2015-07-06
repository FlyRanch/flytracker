function [xx,yy,zz] = xformq_surfMAT(x,y,z,T,q,scale)

%Perform coordinate transformation for a set of points that define a
%surface
if nargin == 4
    q = [0 0 0 1];
    scale = 1;
elseif nargin == 5
    scale = 1;
end

[nrow,ncol] = size(x);

%Make Row;
%T = T(:)';

%Add zeros to work with quaternion math
%pts = [reshape(x,[],1) reshape(y,[],1) reshape(z,[],1) zeros(nrow*ncol,1)];
pts = [x(:) y(:) z(:) zeros(nrow*ncol,1)];

% quaternion math uses formula from Chap 4 of "Intro to theoretical 
% kinematics", by J. McCarthy. Quaternion is in the form [q1 q2 q3 q4] 
% where q4 is the scalar component.

Pplus = [q(4) -q(3) q(2) q(1);
    q(3) q(4) -q(1) q(2);
    -q(2) q(1) q(4) q(3);
    -q(1) -q(2) -q(3) q(4)];

% Pminus = [q(4) q(3) -q(2) -q(1);
%     -q(3) q(4) q(1) -q(2);
%     q(2) -q(1) q(4) -q(3);
%     q(1) q(2) q(3) q(4)];

% Pminusconj = [q(4) -q(3) q(2) -q(1);
%     q(3) q(4) -q(1) -q(2);
%     -q(2) q(1) q(4) -q(3);
%     q(1) q(2) q(3) q(4)];

Pminusconj = Pplus.*[1 1 1 -1 ; 1 1 1 -1 ; 1 1 1 -1; -1 -1 -1 1];
Pplus*Pminusconj

% (Pplus*Pminusconj * pts')' = pts*(Pplus*Pminusconj)' matrix transpose
pts = pts*((Pplus*Pminusconj)');
pts = scale.*pts(:,1:3);

%Add translation
for j=1:3
    pts(:,j) = pts(:,j) + T(j);
end

xx = reshape(pts(:,1),nrow,[]);
yy = reshape(pts(:,2),nrow,[]);
zz = reshape(pts(:,3),nrow,[]);