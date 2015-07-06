function Y = RectVertices(l,w,x0,y0)
%
% Y = RectVertices(l,w,x0,y0)
% Returns the vertices of a rectangle with length 'l', width, 'w' where the
% center in located at (x0,y0)

Y = 0.5.*[-l -w
    -l w
    l w
    l -w
    -l -w];

Y(:,1) = Y(:,1) + x0;
Y(:,2) = Y(:,2) + y0;