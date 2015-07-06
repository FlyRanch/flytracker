% Find the median axis of fly body boundary points
load flygenmod
debug = 0;

c = 4;
npts = size(flygenmod.bodybndy,1);
knotv = oknotP(npts,c,1,0);
T = linspace(-0.5,0.5,300);
N = dbasisP(c,T,npts,knotv);

boundary = N*flygenmod.bodybndy;

extendedboundary = [];
splitBoundaryFly1 = [];
splitBoundaryFly2 = [];
smoothboundary = [];

%order of butterworth filter
filter_order = 4;
%filter parameters definieeren
[b a] = butter( filter_order, 2*1/length(boundary)*50);
%spiegelen van de randen beginpunt is ook eindpunt en dat moet er netjes in
% mirrorboundary=1*boundary;
% %volgorde omdraaien
% omdraaimirrorboundary=mirrorboundary(length(mirrorboundary):-1:1,:);
extendedboundary(1:length(boundary),:)=boundary(1:length(boundary),:);
extendedboundary(length(boundary)+1:2*length(boundary),:)=boundary(1:length(boundary),:);
extendedboundary(2*length(boundary)+1:3*length(boundary),:)=boundary(1:length(boundary),:);
if debug==1
    figure,plot(extendedboundary)
end

LONGsmoothboundary=filtfilt( b, a, extendedboundary);
smoothboundary(1:length(boundary),:)=LONGsmoothboundary(length(boundary)+1:2*length(boundary),:);
%close boundary
% smoothboundary(length(boundary)+1,:)=smoothboundary(1,:);
% figure,plot(smoothboundary)

if debug==1
    figure,imshow(label2rgb(L, @jet, [.5 .5 .5]))
    hold on, plot(smoothboundary(1:length(smoothboundary),2),smoothboundary(1:length(smoothboundary),1), 'w', 'LineWidth', 2);
    hold off
end

figure(1); clf;
plot(smoothboundary(:,1),smoothboundary(:,2),'.-');
axis equal

%start interface clickable start and end point-----------------------------

%click start and end point midline
%find head with mouse
title('zoom head | klick enter in matlab');
zoom on
pause
zoom off
title('click head | enter');
[xStartEnd(1),yStartEnd(1)] = ginput(1);
zoom out
clc;

%find tail with mouse
title('zoom tail | klick enter in matlab');
zoom on
pause
zoom off
title('click tip tail | enter');
[xStartEnd(2),yStartEnd(2)] = ginput(1);
clc;
%end interface clickable start and end point-------------------------------

%I need to determine which point on the smoothboundary is closest
%to the clicked point.  I will use a kdtree instead of a
%'for' loop.
iMinStart = kdtreeidx(smoothboundary,[xStartEnd(1) yStartEnd(1)]);
iMinEnd = kdtreeidx(smoothboundary,[xStartEnd(2) yStartEnd(2)]);

%start determine both outlines Fly devided by snout and tail points-------
if iMinStart<iMinEnd
    splitBoundaryFly1(1:(iMinEnd-iMinStart)+1,:) = smoothboundary(iMinStart:iMinEnd,:);
    splitBoundaryFly2(1:iMinStart,:)=smoothboundary(iMinStart:-1:1,:);
    splitBoundaryFly2(iMinStart+1:iMinStart+1+ ...
        length(smoothboundary)-iMinEnd,:)=smoothboundary(length(smoothboundary):-1:iMinEnd,:);
    if debug==1
        figure,plot(splitBoundaryFly1(:,2),splitBoundaryFly1(:,1),'b');
        hold on
        plot(splitBoundaryFly2(:,2),splitBoundaryFly2(:,1),'r');
        hold off
    end
else
    splitBoundaryFly2(1:(iMinStart-iMinEnd)+1,:)=smoothboundary(iMinStart:-1:iMinEnd,:);
    splitBoundaryFly1(1:length(smoothboundary)-iMinStart+1,:)=smoothboundary(iMinStart:end,:);
    splitBoundaryFly1(length(smoothboundary)-iMinStart+2:...
        length(smoothboundary)-iMinStart+1+iMinEnd,:)=smoothboundary(1:iMinEnd,:);
    if debug==1
        figure,plot(splitBoundaryFly1(:,2),splitBoundaryFly1(:,1),'b');
        hold on
        plot(splitBoundaryFly2(:,2),splitBoundaryFly2(:,1),'r');
        hold off
    end
end
%end determine both outlines Fly devided by snout and tail points---------


% I will take these two sets of pixel boundary locations and
% fit a cubic b-spline to each one that is parameterized over
% the same domain.  The b-spline fit will be my X,Y points.

%Flip the arrays to make them go from tail to snout
splitBoundaryFly1 = flipud(splitBoundaryFly1);
splitBoundaryFly2 = flipud(splitBoundaryFly2);

radpts1 = splitBoundaryFly1;
radpts2 = splitBoundaryFly2;

figure(1); clf
plot(radpts1(:,1),radpts1(:,2),'r.-',radpts2(:,1),radpts2(:,2),'b.-');
hold on;
axis equal
ax1 = gca;

zoom on
pause
zoom off

xrange = get(ax1,'xlim');
yrange = get(ax1,'ylim');


title(['Click along center of fly.  Start at the ' ...
    'tail and move towards head.  Click outside image when done']);

[x,y] = ginput(1);

%Initialize
X = [];
while x > xrange(1) && x < xrange(2) && y > yrange(1) && y < yrange(2)
    X = [X;x y];

    %%%%%%%%%%%%%%%%%%%
    % To approximate length of worm by clicking a bunch of times
    %%%%%%%%%%%%%%%%%%%
    if size(X,1) == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use To plot points on Image
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(X(:,1),X(:,2),'o');
    elseif size(X,1) >= 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use To plot points on Image
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(X(:,1),X(:,2),'o-');
    end
    [x,y] = ginput(1);
end

ctrpts1 = X;

% Err = 1;
% while Err > 1e-2
%     DD = zeros(size(radpts1,1),size(radpts2,1));
%     for i = 1:size(radpts1,1)
%         DD(i,:) = sum((repmat(radpts1(i,:),size(radpts2,1),1) - radpts2).^2,2)';
%     end
% 
%     [dist1,nuidx] = min(DD,[],2);
%     [dist2,nnuidx] = min(DD,[],1);
% 
%     %Get rid of points that have converged ontop of each other
%     ctrpts1 = (radpts1 + radpts2(nuidx,:))./2;
%     dd = sqrt(sum(diff(ctrpts1,1,1).^2,2));
%     %ctrpts1(dd < 0.5,:) = [];
% 
%     ctrpts2 = (radpts1(nnuidx,:) + radpts2)./2;
%     dd = sqrt(sum(diff(ctrpts2,1,1).^2,2));
%     %ctrpts2(dd < 0.5,:) = [];
% 
%     figure(2);
%     hold on;
%     plot(radpts1(:,1),radpts1(:,2),'r.-',radpts2(:,1),radpts2(:,2),'b.-')
%     axis equal;
%     
%     for k = 1:size(radpts1,1)
%         plot([radpts1(k,1) radpts2(nuidx(k),1)],[radpts1(k,2) radpts2(nuidx(k),2)],'g-');
%     end
% 
%     for k = 1:size(radpts2,1)
%         plot([radpts1(nnuidx(k),1) radpts2(k,1)],[radpts1(nnuidx(k),2) radpts2(k,2)],'k-');
%     end
% 
%     hold on; plot(ctrpts1(:,1),ctrpts1(:,2),'k--')
%     hold on; plot(ctrpts2(:,1),ctrpts2(:,2),'g--')
% 
%     Err = sqrt(mean([dist1 ;dist2']));
% 
%     radpts1 = ctrpts1;
%     radpts2 = ctrpts2;
% end
% keyboard
dd = sqrt(sum(diff(ctrpts1,1,1).^2,2));
D = cumsum(dd);
Dist = sum(dd);
T = [0 ; D./Dist ]' - 1;

%If there are large jumps in the points, add a few that are
%linear interpolated
ctrpts1xtmp = interp1(T,ctrpts1(:,1),linspace(-1,0,size(ctrpts1,1)));
ctrpts1ytmp = interp1(T,ctrpts1(:,2),linspace(-1,0,size(ctrpts1,1)));

%Reassign the smoothed points
ctrpts1 = [ctrpts1xtmp' ctrpts1ytmp'];

% Fit spline to centerline
Ytemp = ctrpts1;
npts = round(size(Ytemp,1)/6);
c = 4; %cubic spline

knotv = oknotP(npts,c,1,0);
T = linspace(-0.5,0.5,size(Ytemp,1));
[N,D1,D2] = dbasisP(c,T,npts,knotv);
BB = pinv(N)*Ytemp;
ctrpts = N*BB;

%I will get an initial estimate of the width of the Fly by
%reparameterizing the radius curve
radpts1 = splitBoundaryFly1;
radpts2 = splitBoundaryFly2;

rad1x = interp1(1:size(radpts1,1),radpts1(:,1),linspace(1,size(radpts1,1),size(ctrpts,1)));
rad1y = interp1(1:size(radpts1,1),radpts1(:,2),linspace(1,size(radpts1,1),size(ctrpts,1)));
R = sqrt(sum(([rad1x' rad1y'] - ctrpts).^2,2));

%only positive to account for normal vector facing in the
%wrong direction
R = abs(R);

d = sqrt(sum(diff(ctrpts,1,1).^2,2));
D = cumsum(d);
Dist = sum(d);
s = D;

% s is now my approximate arc length sample parameter.
s = [0 ; s]'; % include 0, and leave out s(end)
len = Dist;
s = s - len/2;

% npts = 15; 
npts = 20; %for body 
knotv = oknot(npts,c,len,0);
N = dbasis(c,s,npts,knotv);
SplinePts = pinv(N)*R;
p0 = mean(N*SplinePts)*ones(size(SplinePts));

%Also, reparameterize the centerline with this parameterization
npts = round(size(Ytemp,1)/6);
knotv = oknotP(npts,c,len,0);
N = dbasisP(c,s,npts,knotv);
BB = pinv(N)*Ytemp;
flygenmod.bodyctr = BB;

PAR.dt =1/6000;
PAR.etamax = 0; 
PAR.c = 4;
PAR.pixpermm = 1;
PAR.numworms = 1;
%Now, refine this estimate by using a nonlinear solver
SplinePts = findFlyWidth(BB,p0,len,[radpts1;radpts2],PAR);
flygenmod.bodyrad = SplinePts;
