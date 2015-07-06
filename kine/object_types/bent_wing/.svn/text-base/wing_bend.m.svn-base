% function [x_out, y_out, z_out] = wing_bend(model,d,a,t)
%
% Given the closed polygon forming the fly wing-shape, wing_bend outputs
% the coordinates (x_out, y_out, z_out) of the same polygon bent theta (t) 
% degrees around an axis defined by d and alpha (a):
%
%       model = 3d coordinates of the model to bend
%       d = percent along major axis, ranges from 0 to 1
%       a = angle bend axis makes with the major axis (defnied ccw from
%           major axis, which runs from wing hinge to tip)
%       t = angle to bend the wing (moving the wing tip end ccw around)
%
%       Angles should be entered in degrees.

function [wing_bent, varargout] = wing_bend(model,d,a,t)

%--------------------------------------------------------------------------
ax_length = 1;
ax_pt_num = 30;
%--------------------------------------------------------------------------

% Convert alpha, theta to radian from polar
a = a*pi/180;
t = t*pi/180;

% move the wing so ax_pt is at origin
model_trans = moveobj3(model,[d 0 0]);    
ind_planar = find(model_trans(3,:) == 0);
ind_not_planar = find(model_trans(3,:) ~= 0); % to add to tip ind

r = linspace(-ax_length,ax_length,ax_pt_num);     
th = a * ones(1,ax_pt_num);             
% ax = polar(t,r);     this would be how to draw the axis
% put the axis in cartesian coords
[axis_x, axis_y] = pol2cart(th,r);

% To find points to move, put all points in polar coordinates using new origin
[model_t, model_r] = cart2pol(model_trans(1,ind_planar), model_trans(2,ind_planar));

% Get all the points with theta from -(180-a) to a; ind is the vector of indices of points to move
ind = find(model_t > (a -pi) & model_t < a);
ind = union(ind_planar(ind), ind_not_planar);
tip = model_trans(:,ind);

% % find the tip index within the tip section
% rel_tip_ind = find(model(1,:) == 1 & model(2,:) == 0 & model(3,:) == 0);
% if isempty(ind), rel_tip_ind = rel_tip_ind;
% else rel_tip_ind = rel_tip_ind - ind(1) + 1;
% end
% !!!!!!!!! For now assume that each bend adds 1 to the tip index

% Find intersection of axis line and wing perimeter
% Because now bending wing twice, need to only find intersection with
% points that are still planar, i.e. z = 0
[int_x, int_y] = polyxpoly(model_trans(1,ind_planar), model_trans(2,ind_planar), axis_x, axis_y);

% Bend the wing
tip_rot = rot3D(tip,a,0,0);             % rotate wing about z-axis to align bending axis with x-axis
tip_rot_bend = rot3D(tip_rot,0,0,t);    % bend model in 3d around x-axis
tip_bend = rot3D(tip_rot_bend,-a,0,0);  % unrotate around z-axis

% unmove from new origin
tip = moveobj3(tip_bend,-[d 0 0]);      

% move axis back
[axis_x axis_y] = move(axis_x,axis_y,-[d 0]);

% move intersection points back
[int_x int_y] = move(int_x,int_y,-[d 0]);
int_z = zeros(size(int_x));

% Get total coordinates for bent wing
if length(ind) == length(model)
    x = tip(1,:);
    y = tip(2,:);
    z = tip(3,:);
    tip_ind = rel_tip_ind;
    hinge_ind = find(tip(1,:) == 0 & tip(2,:) == 0 & tip(3,:) == 0);
    int_ind = hinge_ind;
elseif isempty(ind)
    x = model(1,:);
    y = model(2,:);
    z = zeros(size(model(1,:)));
    tip_ind = find(model(1,:) == 1 & model(2,:) == 0 & model(3,:) == 0);
    hinge_ind = find(model(1,:) == 0 & model(2,:) == 0 & model(3,:) == 0);
    int_ind = tip_ind;
else
    x = model(1,1:ind(1)-1);                                            % first segment is points from hinge to just before bend
    y = model(2,1:ind(1)-1);
    z = model(3,1:ind(1)-1);
    
    dist_1 = pt_dist3(x(end),y(end),z(end),int_x(1),int_y(1),int_z(1)); % find the intersection point closest to last point of first section
    dist_2 = pt_dist3(x(end),y(end),z(end),int_x(2),int_y(2),int_z(2));
    
    if dist_1 < dist_2                                                  % closest intersection is the next point
        x(end+1) = int_x(1);
        y(end+1) = int_y(1);
        z(end+1) = int_z(1);
        other = 2;
    elseif dist_2 < dist_1
        x(end+1) = int_x(2);
        y(end+1) = int_y(2);
        z(end+1) = int_z(2);
        other = 1;
    else
        warndlg('Correct intersection point cannot be found')
        return
    end

    int_ind(1) = length(x);                                             % index of one intersection point
    
    x = [x tip(1,:)];                                                   % next segment is the wing tip
    y = [y tip(2,:)];
    z = [z tip(3,:)];

    %     tip_ind = int_ind(1) + rel_tip_ind;                           % index of tip point
        
    x(end+1) = int_x(other);                                            % add other intersection point
    y(end+1) = int_y(other);
    z(end+1) = int_z(other);

    int_ind(2) = length(x);                                             % index of other intersection point

    x = [x model(1,(ind(end)+1):end)];                                  % add on other wing base section
    y = [y model(2,(ind(end)+1):end)];
    z = [z model(3,(ind(end)+1):end)];

    %     hinge_ind = find(model(1,:) == 0 & model(2,:) == 0 & model(3,:) == 0);

end

wing_bent = [x;y;z];
varargout{1} = int_ind;

% varargout{1} = hinge_ind;
% varargout{2} = tip_ind;
% varargout{3} = int_ind;

%--------------------------------------------------------------------------
function d = pt_dist3(x1,y1,z1,x2,y2,z2)

d = sqrt( ( x2 - x1 ).^2 + ( y2 - y1 ).^2 + ( z2 - z1 ).^2 );

%--------------------------------------------------------------------------
function [x, y] = move(x_old,y_old,new_origin)

x = x_old - new_origin(1);
y = y_old - new_origin(2);

%--------------------------------------------------------------------------
function [obj_coords_new] = moveobj3(obj_coords,new_origin)

for i = 1:size(obj_coords,1)
    obj_coords_new(i,:) = obj_coords(i,:) - new_origin(i);
end