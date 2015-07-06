% Function [wing_bent, hinge, tip, int, len] = WING_DOUBLE_BEND(wing, th1, th2)
% 
% CALLING FUNCTION: obj_program for bent_wing_2 types
% ACTIONS: Uses wing_bend to bend the wing model in two places, then rotates and scales
%          so distance/orientation of vector from wing hinge to bent wing tip is the
%          same as that for unbent wing:
%
% Given the 3D closed polygon forming the fly wing-shape, wing_bend outputs
% the coordinates, wing_bent, of the same polygon bent theta (th1) 
% degrees around an axis defined by d=0.13 and a = 78 degrees, theta (th2) degrees
% around an axis defined by d = 0.19 and a = 37 degrees, and then rotated and
% scaled to align with x-axis.
%
%       wing = 3d coordinates of the model to bend
%       th1 = angle to bend the wing (moving the wing tip end ccw around)
%       th2 = second axis
%       Angles should be entered in degrees.
%
% PARENT PROGRAM: Kine_v2_1
% LAST MODIFIED: October 12, 2005 by gwyneth


% function [wing_bent, hinge, tip, int, len] = wing_double_bend(wing,  th1, th2)
function [wing_bent, ax1, ax2] = wing_double_bend(wing,  th1, th2)

%--------------------------------------------------------------------------
% HARD-WIRED axis definitions
d1 = 0.13;
a1 = 90;%78;

d2 = 0.21;%0.19;
a2 = 48;%37;

% HARD-WIRED TIP AND HINGE INDICES
hinge = 1;
tip = 21;
%--------------------------------------------------------------------------

% figure, hold on
% plot3(wing(1,:),wing(2,:),wing(3,:))

% Bend frame around second axis (farther from hinge)
[wing_bent, ax1] = wing_bend( wing, d2, a2, th2 );
tip = tip + 1;

% plot3(wing_bent(1,:),wing_bent(2,:),wing_bent(3,:),'m')

% Bend frame around first axis (closer to hinge)
[wing_bent, ax2] = wing_bend( wing_bent, d1, a1, th1 );
tip = tip + 1;

% plot3(wing_bent(1,:),wing_bent(2,:),wing_bent(3,:),'r')

% Rotate bent frame to align hinge-tip vectors
[th, phi, r] = cart2sph(wing_bent(1,tip),wing_bent(2,tip),wing_bent(3,tip));

wing_bent = rot3d(wing_bent,th,0,0);
wing_bent = rot3d(wing_bent,0,phi,0);

% Scale (important to do this before next rotation)
f = 1/r;
wing_bent = f*wing_bent;

% % Rotate about x axis so bending axis is parallel to xy-plane
% a = a2*pi/180;
% if wing_bent(3,int(2)) - wing_bent(3,int(1)) > 10^-10
%     dy = wing_bent(2,int(2)) - wing_bent(2,int(1));
%     dz = wing_bent(3,int(2)) - wing_bent(3,int(1));
%     [th, r] = cart2pol(dy,dz);
%     
%     if a >= pi/2, th = th+pi; end
%     if t <= 0, th = th+pi; end
%     
% elseif wing_bent(3,int(1)) - wing_bent(3,int(2)) > 10^-10
%     dy = wing_bent(2,int(1)) - wing_bent(2,int(2));
%     dz = wing_bent(3,int(1)) - wing_bent(3,int(2)); 
%     [th, r] = cart2pol(dy,dz);
%     
%     if a >= pi/2, th = th+pi; end
%     if t <= 0, th = th+pi; end
%     
% else
%     if wing_bent(3,int(1)) < 0  
%         if t < 0, th = 0;
%         elseif t > 0, th = pi;
%         end
%     elseif wing_bent(3,int(1)) > 0 
%         if t < 0, th = pi;
%         elseif t > 0, th = 0;
%         end
%     else 
%         th = 0;
%     end
% end
% 
% wing_bent = rot3d(wing_bent,0,0,-th);

% % Find bent wing length
% pt = [d; 0; wing_bent(3,int(1))];
% 
% len1 = pt_dist3(wing_bent(1,hinge(1)),wing_bent(2,hinge(1)),wing_bent(3,hinge(1)),pt(1),pt(2),pt(3));
% len2 = pt_dist3(wing_bent(1,tip),wing_bent(2,tip),wing_bent(3,tip),pt(1),pt(2),pt(3));
% 
% len = len1 + len2;