function [mod_rot] = roteulMstar(model_3d,q);

[phi, th, psi] = quat2eul_Mstar(q);

phi*180/pi

model_3d = [0 0 1; 0 1 0; -1 0 0]*model_3d; % starting position on negative z-axis

% m_st = [0 0 1; 0 1 0; -1 0 0];
% m_phi = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
% m_th = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
% m_psi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
% 
% v1 = m_st*model_3d;
% v2 = m_phi*v1;
% v3 = m_th*v2;
% v4 = m_psi*v3;
% 
% mod_rot = v4;

%pi/2 rotation
M_st = [ -cos(psi)*sin(th)*cos(phi)-sin(psi)*sin(phi)   cos(psi)*sin(th)*sin(phi)-sin(psi)*cos(phi) cos(psi)*cos(th);...
         -sin(psi)*sin(th)*cos(phi)+cos(psi)*sin(phi)   sin(psi)*sin(th)*sin(phi)+cos(psi)*cos(phi) sin(psi)*cos(th);...
         -cos(th)*cos(phi)                              cos(th)*sin(phi)                            -sin(th)        ];

% %-pi/2 rotation
% M_st = [ cos(psi)*sin(th)*cos(phi)+sin(psi)*sin(phi)   cos(psi)*sin(th)*sin(phi)-sin(psi)*cos(phi) -cos(psi)*cos(th);...
%          sin(psi)*sin(th)*cos(phi)-cos(psi)*sin(phi)   sin(psi)*sin(th)*sin(phi)-cos(psi)*cos(phi) -sin(psi)*cos(th);...
%          cos(th)*cos(phi)                              cos(th)*sin(phi)                            sin(th)        ];
     
mod_rot = M_st*model_3d;


% figure
% plot3(model_3d(1,:),model_3d(2,:),model_3d(3,:),'b')
% hold on
% plot3(mod_rot(1,:),mod_rot(2,:),mod_rot(3,:),'m')