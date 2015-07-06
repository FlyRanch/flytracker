function [x, y] = move(x_old,y_old,new_origin)

x = x_old - new_origin(1);
y = y_old - new_origin(2);