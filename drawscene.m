function drawscene(X,C,R,fig,ctypehint,scenetitle,camsId)

% drawscene ... plots calibration points, cameras and their viewing axes
%
% drawscene(X,C,R,fig,ctypehint,scenetitle)
%
% X ............ 4xPOINTS matrix containg POINTS object points
% C ............ 3xCAMS matrix containing the camera centers (in world coord.)
% R ............ 3*CAMSx3 matrix containing camera rotation matrices
%                (needed for drawing the viewing axes)
% fig .......... figure handle (defaults to 1)
% ctypehint .... calibration object type of X (defaults to 'cloud')
% scenetitle ... title of the plot (defaults to '')
% camsIs ....... 1xCAMS vector with cameas Id (default is 1:CAMS

% $Author: svoboda $
% $Revision: 2.0 $
% $Id: drawscene.m,v 2.0 2003/06/19 12:07:03 svoboda Exp $
% $State: Exp $

POINTS = size(X,2);
CAMS   = size(C,2);

if nargin < 7
  camsId = [1:CAMS];
end

if (nargin < 3)
  error('not enough input arguments');
end
if (nargin < 5)
  scenetitle = '';
end
if (nargin < 4)
  ctypehint = 'cloud';
end

figure(fig); clf
title(scenetitle)
grid on
axis equal

% plot camera positions (blue)
drawcloud(C,fig,'b');

% plot calibration object (red)
if ~strcmp(ctypehint,'none')
    drawobject(X,ctypehint,fig,'r');
end
    drawobject(X,ctypehint,fig,'r');
% Mean of all points
centroid = mean(X(1:3,:)');

%Plot world reference frame
ref = 20*eye(3);
org = zeros(3,3);
%quiver3(org(:,1),org(:,2),org(:,3),ref(:,1),ref(:,2),ref(:,3));
%hold on;

% plot viewing axes
for i=1:CAMS
    camframe = [R(3*i-2:3*i,:)' C(1:3,i)] * [org ref;ones(1,6)];
    xax = [camframe(1:3,1) camframe(1:3,4)]'; 
    yax = [camframe(1:3,2) camframe(1:3,5)]';
    zax = [camframe(1:3,3) camframe(1:3,6)]';
    %plot3(xax(:,1),xax(:,2),xax(:,3),'r.-',yax(:,1),yax(:,2),yax(:,3),'g.-',...
    %zax(:,1),zax(:,2),zax(:,3),'b.-');
    plot3(C(1,i),C(2,i),C(3,i),'o','markersize',15,'markerfacecolor','k')
    hold on;
    
  axis_dir = R(3*i,:); % 3rd row of i-th rotation matrix
%   axis_dir = R(3*i-2:3*i,3)'; % 3rd row of i-th rotation matrix
  axis_len = 0.6*norm(C(1:3,i)-centroid');  
  endpoint = C(1:3,i)+axis_len*axis_dir';
  %line([C(1,i),endpoint(1)],[C(2,i),endpoint(2)],[C(3,i),endpoint(3)]);
  text(C(1,i),C(2,i),C(3,i),sprintf('%4d',camsId(i)),'Color','k');
end


