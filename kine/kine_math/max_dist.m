% function [ ind_1, ind_2 ] = max_dist(x,y,z)
%
% Given a set of 3d coordinates, finds the two points that are fartheset
% apart.  The output is the two index numbers of these points in the x, y,
% z matrices.

function [ind_1, ind_2] = max_dist(x,y,z)

for i = 1:length(x)

    dist_table(i,:) = pt_dist3(x(i),y(i),z(i),x,y,z);   % Distance from the ith point to all the other points
    
end

[ max1, ind_1 ] = max(max(dist_table));

[ max2, ind_2 ] = max(dist_table(ind_1,:));