function [mod_rot] = roteulMstar8(model_3d,q)

[eul] = quat2eul_rotquat_Mstar5(q); % output in DEGREES!!

phi = eul(1)*pi/180;
th = eul(2)*pi/180;
psi = eul(3)*pi/180;

M_st = [ cos(psi)*cos(th) -cos(phi)*sin(psi)+sin(phi)*sin(th)*cos(psi)   sin(phi)*sin(psi)+cos(phi)*sin(th)*cos(psi) ;...
           sin(psi)*cos(th)  cos(phi)*cos(psi)+sin(phi)*sin(th)*sin(psi)  -sin(phi)*cos(psi)+cos(phi)*sin(th)*sin(psi) ;...
           -sin(th)          cos(th)*sin(phi)                              cos(th)*cos(phi)                            ];


mod_rot = M_st*model_3d;

