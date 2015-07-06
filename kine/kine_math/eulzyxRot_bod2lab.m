% Function EULZYXROT_BOD2LAB(phi, theta, psi, vb)
% 
% CALLING FUNCTION: obj_function (frame)
% ACTIONS: Rotate the vector vb "forward" using Euler angles
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: September 26, 2007 by gwyneth

function vf = eulzyxRot_bod2lab(phi, theta, psi, vb)

ph = phi;
th = theta;
ps = psi;

Mx = [ [1 0 0]; [0 cos(ph) -sin(ph)]; [0 sin(ph) cos(ph)] ];
My = [ [cos(th) 0 sin(th)]; [0 1 0]; [-sin(th) 0 cos(th)] ];
Mz = [ [cos(ps) -sin(ps) 0]; [sin(ps) cos(ps) 0]; [0 0 1] ];

vf = Mz*My*Mx*vb;