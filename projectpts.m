function Y = projectpts(X,cam)
%
%Y = projectpts(X,cam)
%Calculate the pixel locations of the points in X given the camera model in
%cam.

%Let's transform the world coordinates into the camera frame
campts = [cam.R -cam.R*cam.T]*[X ones(size(X,1),1)]';

%Normalize the point coordinates by perspective.
xn = campts(1:2,:)./[campts(3,:);campts(3,:)];


r2 = sum(xn.^2,1);
        
k1 = cam.k;

KK = [cam.f(1) 0 cam.u0(1)
    0 cam.f(2) cam.u0(2)
    0 0 1];

xd = xn.* repmat(1 + k1*r2,2,1);

Y = KK*[xd ; ones(1,size(xd,2))];
Y = Y(1:2,:)';