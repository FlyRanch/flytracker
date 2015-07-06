% Called by rot3D

function Q = q2Q(q)

e_0 = q(1);
e_x = q(2);
e_y = q(3);
e_z = q(4);

Q = [ e_0^2+e_x^2-e_y^2-e_z^2 2*e_x*e_y-2*e_0*e_z     2*e_x*e_z+2*e_0*e_y;...
      2*e_x*e_y+2*e_0*e_z     e_0^2-e_x^2+e_y^2-e_z^2 2*e_y*e_z-2*e_0*e_x;...
      2*e_x*e_z-2*e_0*e_y     2*e_y*e_z+2*e_0*e_x     e_0^2-e_x^2-e_y^2+e_z^2 ];