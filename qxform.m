function RR = qxform(q,x)

% RR = qxform(q,x) transforms the point x with quaternion q.  Uses
% formula from "Intro to theoretical kinematics", by J. McCarthy.
% Quaternion is in the form [q1 q2 q3 q0] where q0 is the scalar
% component.

%Rearrange q to work with this construction
q = q(:);


Pplus = [q(4) -q(3) q(2) q(1);
	 q(3) q(4) -q(1) q(2);
	 -q(2) q(1) q(4) q(3);
	 -q(1) -q(2) -q(3) q(4)];

Pminusconj = [q(4) -q(3) q(2) -q(1);
	 q(3) q(4) -q(1) -q(2);
	 -q(2) q(1) q(4) -q(3);
	 q(1) q(2) q(3) q(4)];

RR = Pplus*Pminusconj*[x;0];

RR = RR(1:3,:);
