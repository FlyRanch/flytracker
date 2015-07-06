% Function [wing_bent, hinge, tip, int, len] = WING_BEND2(wing, d, a, t)
% 
% CALLING FUNCTION: obj_program for bent_wing types
% ACTIONS: Uses wing_bend to bend the wing model, then rotates and scales
%          so distance/orientation of vector from wing hinge to bent wing tip is the
%          same as that for unbent wing:
%
% Given the 3D closed polygon forming the fly wing-shape, wing_bend outputs
% the coordinates, wing_bent, of the same polygon bent theta (t) 
% degrees around an axis defined by d and alpha (a), and then rotated and
% scaled to align with x-axis.
%
%       wing = 3d coordinates of the model to bend
%       d = percent along major axis, ranges from 0 to 1
%       a = angle bend axis makes with the major axis (defnied ccw from
%           major axis, which runs from wing hinge to tip)
%       t = angle to bend the wing (moving the wing tip end ccw around)
%
%       Angles should be entered in degrees.
%
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: November 16, 2004 by gwyneth


function [wing_bent, hinge, tip, int, len] = wing_bend2(wing, d, a, t)

% Bend frame
[wing_bent, hinge, tip, int] = wing_bend( wing, d, a, t );

% Rotate bent frame to align hinge-tip vectors
a = a*pi/180;
[th, phi, r] = cart2sph(wing_bent(1,tip),wing_bent(2,tip),wing_bent(3,tip));

wing_bent = rot3d(wing_bent,th,0,0);
wing_bent = rot3d(wing_bent,0,phi,0);

% Scale (important to do this before next rotation)
f = 1/r;
wing_bent = f*wing_bent;

% Rotate about x axis so bending axis is parallel to xy-plane
if wing_bent(3,int(2)) - wing_bent(3,int(1)) > 10^-10
    dy = wing_bent(2,int(2)) - wing_bent(2,int(1));
    dz = wing_bent(3,int(2)) - wing_bent(3,int(1));
    [th, r] = cart2pol(dy,dz);
    
    if a >= pi/2, th = th+pi; end
    if t <= 0, th = th+pi; end
    
elseif wing_bent(3,int(1)) - wing_bent(3,int(2)) > 10^-10
    dy = wing_bent(2,int(1)) - wing_bent(2,int(2));
    dz = wing_bent(3,int(1)) - wing_bent(3,int(2)); 
    [th, r] = cart2pol(dy,dz);
    
    if a >= pi/2, th = th+pi; end
    if t <= 0, th = th+pi; end
    
else
    if wing_bent(3,int(1)) < 0  
        if t < 0, th = 0;
        elseif t > 0, th = pi;
        end
    elseif wing_bent(3,int(1)) > 0 
        if t < 0, th = pi;
        elseif t > 0, th = 0;
        end
    else 
        th = 0;
    end
end

wing_bent = rot3d(wing_bent,0,0,-th);

% Find bent wing length
pt = [d; 0; wing_bent(3,int(1))];

len1 = pt_dist3(wing_bent(1,hinge(1)),wing_bent(2,hinge(1)),wing_bent(3,hinge(1)),pt(1),pt(2),pt(3));
len2 = pt_dist3(wing_bent(1,tip),wing_bent(2,tip),wing_bent(3,tip),pt(1),pt(2),pt(3));

len = len1 + len2;