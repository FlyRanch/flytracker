function [wing_bent, hinge, tip, int, len] = testbend(d, a, t);

load('C:\MATLAB7\work\kine_v2_1\models\wing_fly_bent.mat')

wing = [coords; zeros(1,length(coords))];

% d = 0.1;
% a = 100;
% t = -30;

labels = anchor_array(:,1);           % get from wingbend
tip_ind = strmatch('tip',labels);
tip = anchor_array{tip_ind,2};

figure
view(-30,20)
set(gca, 'XLim', [-0.4 1.4])
set(gca, 'YLim', [-0.4 0.4])
set(gca, 'ZLim', [-0.4 0.4])

hold on
plot3(wing(1,:), wing(2,:), wing(3,:))
plot3(wing(1,tip), wing(2,tip), wing(3,tip),'ro')
plot3(wing(1,1), wing(2,1), wing(3,1),'ro')

pause

[wing_bent, hinge, tip, int] = wing_bend(wing,d,a,t);
plot3(wing_bent(1,:),wing_bent(2,:),wing_bent(3,:),'c:')
%tip = tip + 1; %%%FIX!!!!!!! b

pause

% cart2pol in the xz plane (do phi, theta, phi to get rid of imbiguous
% thetas)
a = a*pi/180;
[th, phi, r] = cart2sph(wing_bent(1,tip),wing_bent(2,tip),wing_bent(3,tip));
% th*180/pi
% phi*180/pi

wing_bent = rot3D(wing_bent,th,0,0);
wing_bent = rot3D(wing_bent,0,phi,0);

plot3(wing_bent(1,:),wing_bent(2,:),wing_bent(3,:),'g:')
%scale
f = 1/r;
wing_bent = f*wing_bent;
%plot3(wing_bent(1,:),wing_bent(2,:),wing_bent(3,:),'r:')
%pause

% now need to rotate about x axis, take random point from either side of
% tip, make elevation of vector between these two points zero

len = 1;

plot3(wing_bent(1,int(1)), wing_bent(2,int(1)), wing_bent(3,int(1)),'go')
plot3(wing_bent(1,int(2)), wing_bent(2,int(2)), wing_bent(3,int(2)),'go')
pause

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

% th*180/pi
wing_bent = rot3d(wing_bent,0,0,-th);

plot3(wing_bent(1,:),wing_bent(2,:),wing_bent(3,:),'r')
plot3(wing_bent(1,tip), wing_bent(2,tip), wing_bent(3,tip),'mo')
plot3(wing_bent(1,1), wing_bent(2,1), wing_bent(3,1),'mo')

