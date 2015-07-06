function [x, y, z] = move3(x_old,y_old,z_old,new_origin)

x = x_old - new_origin(1);
y = y_old - new_origin(2);
z = z_old - new_origin(3);