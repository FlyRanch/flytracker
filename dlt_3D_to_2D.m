function [ u ] = dlt_3D_to_2D( DLT,threeDpts)
% DLT_3D_TO_2D: use DLT transformation, given by 11 coefficients in the
% array DLT, to calculate the point in the image plane (u,v) corresponding 
% to the 3d point (x,y,z).

x = threeDpts(:,1);
y = threeDpts(:,2);
z = threeDpts(:,3);

denom = DLT(9)*x + DLT(10)*y + DLT(11)*z + 1;

u(:,1) = ( DLT(1)*x + DLT(2)*y + DLT(3)*z + DLT(4) )./denom;

u(:,2) = ( DLT(5)*x + DLT(6)*y + DLT(7)*z + DLT(8) )./denom;


