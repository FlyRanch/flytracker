%Rotate points [[x1;y1;z1],[x2;y2;z2],...] by angle theta around an arbitrary axis r. Return the rotated points (again each point is a column).

%   Positive angles are anticlockwise looking down the axis
%   towards the origin.

function q=arbitraryRotate(p,theta,r)
%q=[0 0 0];
q=zeros(size(p));
r=r/len(r);
costheta = cos(theta);
sintheta = sin(theta);

q(1,:) = q(1,:) + (costheta + (1 - costheta) * r(1) * r(1)) * p(1,:);
q(1,:) = q(1,:) + ((1 - costheta) * r(1) * r(2) - r(3) * sintheta) * p(2,:);
q(1,:) = q(1,:) + ((1 - costheta) * r(1) * r(3) + r(2) * sintheta) * p(3,:);

q(2,:) = q(2,:) + ((1 - costheta) * r(1) * r(2) + r(3) * sintheta) * p(1,:);
q(2,:) = q(2,:) + (costheta + (1 - costheta) * r(2) * r(2)) * p(2,:);
q(2,:) = q(2,:) + ((1 - costheta) * r(2) * r(3) - r(1) * sintheta) * p(3,:);

q(3,:) = q(3,:) + ((1 - costheta) * r(1) * r(3) - r(2) * sintheta) * p(1,:);
q(3,:) = q(3,:) + ((1 - costheta) * r(2) * r(3) + r(1) * sintheta) * p(2,:);
q(3,:) = q(3,:) + (costheta + (1 - costheta) * r(3) * r(3)) * p(3,:);

