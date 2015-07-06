% Function EULZYXROT_LAB2BOD(phi, theta, psi, vf)
% 
% CALLING FUNCTION: obj_function (frame)
% ACTIONS: Rotate the vector vf "backward" using Euler angles
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: September 26, 2007 by gwyneth

function vb = eulzyxRot_lab2bod(phi, theta, psi, vf)

ph = phi;
th = theta;
ps = psi;

Mx = [ [1 0 0]; [0 cos(ph) -sin(ph)]; [0 sin(ph) cos(ph)] ];
My = [ [cos(th) 0 sin(th)]; [0 1 0]; [-sin(th) 0 cos(th)] ];
Mz = [ [cos(ps) -sin(ps) 0]; [sin(ps) cos(ps) 0]; [0 0 1] ];

vb = Mx'*My'*Mz'*vf;