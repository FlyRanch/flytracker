function cam = dlt2cam(L)
% cam = dlt2cam(L)
%
% This function calculates the transformation matrix, translation,
% camera scaling, and focal length of camera from the 11 DLT parameters
%
%Equations are obtained from http://www.kwon3d.com/theory/dlt/dlt.html
%

%make column
L = L(:);

C = inv([L(1:3)'
    L(5:7)'
    L(9:11)'])*[-L(4) ; -L(8) ; -1];

D2 = 1/sum(L(9:11).^2);

u0 = D2*(L(1)*L(9) + L(2)*L(10) + L(3)*L(11));
v0 = D2*(L(5)*L(9) + L(6)*L(10) + L(7)*L(11));

du2 = D2*( (u0*L(9)-L(1))^2 + (u0*L(10)-L(2))^2 + (u0*L(11)-L(3))^2 );
dv2 = D2*( (v0*L(9)-L(5))^2 + (v0*L(10)-L(6))^2 + (v0*L(11)-L(7))^2 );

du = sqrt(du2);
dv = sqrt(dv2);

%Calculate rotation matrix
D = sqrt(D2);

R = D.*[ [u0*L(9)-L(1) u0*L(10)-L(2) u0*L(11)-L(3)]./du
    [v0*L(9)-L(5) v0*L(10)-L(6) v0*L(11)-L(7)]./dv
    L(9:11)'];

if det(R) < 0
    R = -R;
end

cam.C = C;
cam.R = R;
cam.K = [-du 0 u0 ; 0 -dv v0 ; 0 0 1];
cam.P = [-du 0 u0 ; 0 -dv v0 ; 0 0 1]*[R -R*C]; 
