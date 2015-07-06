function [mod_rot] = roteulMstar(model_3d,q);

[eul] = quat2eul_Mstar3(q); % output in DEGREES!!
phi = eul(1)*pi/180;
th = eul(2)*pi/180;
psi = eul(3)*pi/180-pi;

M_st = [-sin(th)*cos(psi) -cos(psi)*cos(th)*sin(phi)-sin(psi)*cos(phi) cos(psi)*cos(th)*cos(phi)-sin(psi)*sin(phi);...
        -sin(th)*sin(psi) -sin(psi)*cos(th)*sin(phi)+cos(psi)*cos(phi) sin(psi)*cos(th)*cos(phi)+cos(psi)*sin(phi);...
        -cos(th)                    sin(th)*sin(phi)                           -sin(th)*cos(phi)                 ];

     
mod_rot = M_st*model_3d;

