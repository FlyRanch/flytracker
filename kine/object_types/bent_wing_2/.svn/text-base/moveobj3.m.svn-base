% [obj_coords_new] = moveobj3(obj_coords,new_origin)
%
% obj_coords = coordinates of object to move, each coordinate axis is a row
% new_origin = new origin to move object to

function [obj_coords_new] = moveobj3(obj_coords,new_origin)


for i = 1:size(obj_coords,1)
    obj_coords_new(i,:) = obj_coords(i,:) - new_origin(i);
end