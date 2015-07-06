% Function vb = QUATROT_LAB2BOD(q,vf)
% 
% CALLING FUNCTION: obj_program (frame)
% ACTIONS: Frame rotation (backwards) of the vector vf to body frame, vb
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: September 26, 2007 by gwyneth

function vb = quatRot_lab2bod(q,vf)

e_0 = q(1);
e_x = q(2);
e_y = q(3);
e_z = q(4);

Q = [ e_0^2+e_x^2-e_y^2-e_z^2 2*e_x*e_y-2*e_0*e_z     2*e_x*e_z+2*e_0*e_y;...
      2*e_x*e_y+2*e_0*e_z     e_0^2-e_x^2+e_y^2-e_z^2 2*e_y*e_z-2*e_0*e_x;...
      2*e_x*e_z-2*e_0*e_y     2*e_y*e_z+2*e_0*e_x     e_0^2-e_x^2-e_y^2+e_z^2 ];
  
Q = Q';

vb = Q*vf;