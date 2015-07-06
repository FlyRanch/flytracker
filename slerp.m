function Q = slerp(q1,q2,t)

%Q = slerp(q1,q2,t) performs spherical interpolation between q1 and
%q2.  It returns the interpolated quaternion Q.  t is in [0 1].
%Function accepts vector arguments for t and returns a Nx4 matrix
%of quaternions.
q1 = q1(:);
q2 = q2(:);

B = acos(q1'*q2);

Q = zeros(4,length(t));

for n = 1:length(t)
  Q(:,n) = q1*(sin(B*(1-t(n))) / sin(B)) + ...
	   q2*(sin(B*t(n)) / sin(B));
end

