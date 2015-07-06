function Y = quat2mat(X)
% Y = quat2mat(X) converts rotations in SO(3) between matrix and
% quaternion representation.  If one in input, the other is
% output.  Based on "Animating Rotation with Quaternion Curves", by
% K. Shoemake, 1985.
eps = realmin;
dim = prod(size(X));
switch dim
 case 9 %A rotation matrix is the input
  w2 = (1 + X(1,1) + X(2,2) + X(3,3)) / 4;
  if w2 > eps
    w = sqrt(w2);
    x = (X(2,3) - X(3,2)) / (4*w);
    y = (X(3,1) - X(1,3)) / (4*w);
    z = (X(1,2) - X(2,1)) / (4*w);
  else
    w = 0;
    x2 = -(X(2,2) + X(3,3))/2;
    if x2 > eps
      x = sqrt(x2);
      y = X(1,2)/(2*x);
      z = X(1,3)/(2*x);
    else
      x = 0;
      y2 = (1 - X(3,3))/2;
      if y2 > eps
	y = sqrt(y2);
	z = X(2,3)/(2*y);
      else
	y = 0; 
	z = 1;
      end
    end
  end
  
  Y = [x y z w];
  
 case 4 %A quaternion is the input
  % Check to make sure that it is unit magnitude
%   if norm(X) ~= 1
%     error(['the quaternion does not have unit length.  Length = ' ...
% 	   'num2str(norm(X))']);
%   end
  
  w = X(4);
  x = X(1);
  y = X(2);
  z = X(3);
  
  Y = [1 - 2*y^2 - 2*z^2, 2*x*y+2*w*z, 2*x*z-2*w*y
       2*x*y-2*w*z, 1 - 2*x^2 - 2*z^2, 2*y*z+2*w*x
       2*x*z+2*w*y, 2*y*z-2*w*x, 1 - 2*x^2 - 2*y^2];
end


