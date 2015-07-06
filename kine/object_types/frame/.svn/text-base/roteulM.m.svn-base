function [mod_rot] = roteulMstar(model_3d,q);

[eul] = quat2eul_vout(q);
phi = eul(1)*pi/180;
th = eul(2)*pi/180;
psi = eul(3)*pi/180;

% m_phi = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
% m_th = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
% m_psi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
% 
% v1 = model_3d;
% v2 = m_phi*v1;
% v3 = m_th*v2;
% v4 = m_psi*v3;
% 
% mod_rot = v4;

M = [cos(th)*cos(psi)   sin(phi)*sin(th)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(th)*cos(psi)+sin(phi)*sin(psi);...
     cos(th)*sin(psi)   sin(phi)*sin(th)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(th)*sin(psi)-sin(phi)*cos(psi);...
     -sin(th)           sin(phi)*cos(th)                            cos(phi)*cos(th)                           ];
 
mod_rot = M*model_3d;