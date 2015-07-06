% test rotate


obj = [0 1; 0 0; 0 0];
tip = 2;

figure
hold on

plot3(obj(1,:), obj(2,:), obj(3,:))
pause

obj_rot = rot3D(obj, -30*pi/180, -30*pi/180, 0);
plot3(obj_rot(1,:), obj_rot(2,:), obj_rot(3,:),'c')
pause

[theta, phi, r] = cart2sph(obj_rot(1,tip),obj_rot(2,tip), obj_rot(3,tip));

theta*180/pi
phi*180/pi
r

obj_rot_back = rot3D(obj_rot, theta, 0, 0);
obj_rot_back = rot3D(obj_rot_back, 0, phi, 0);
plot3(obj_rot_back(1,:), obj_rot_back(2,:), obj_rot_back(3,:),'r')