function [mod_rot] = roteulMstar(model_3d,q);

[eul] = quat2eul_Mstar5(q); % output in DEGREES!!
phi = eul(1)*pi/180;
th = eul(2)*pi/180;
psi = eul(3)*pi/180;

%model_3d = [0 0 1; 0 1 0; -1 0 0]*model_3d; % starting position on negative z-axis

M_st = [  cos(psi)*cos(th)*cos(phi)-sin(psi)*sin(phi)  -cos(psi)*cos(th)*sin(phi)-sin(psi)*cos(phi)  sin(th)*cos(psi);...
    sin(psi)*cos(th)*cos(phi)+cos(psi)*sin(phi)  -sin(psi)*cos(th)*sin(phi)+cos(psi)*cos(phi)  sin(th)*sin(psi);...
    -sin(th)*cos(phi)                              sin(th)*sin(phi)                             cos(th)        ];


mod_rot = M_st*model_3d;

