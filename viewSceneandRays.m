% viewSceneandRays.m
%
% This file draws a picture of the experimental setup and projection rays
% calculated from 'feat_detect.m'
close all;clc;clear all

PAR.videopath = '../FlyVideo/';
PAR.filetag = 'exp035';

PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.stub = [PAR.filetag];

fignum = 1;
frame = 134;

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .75;

cptsTag = 0;   %Set to 1 to see the corresponding points

PAR.pixpermm = 1;
PAR.numfly = 1;
%Number of parameters of the model (i.e. 8 control points)
PAR.mdlpar = 15*ones(1,PAR.numfly);
PAR.modelfun_H = @modcurvesplineP;
PAR.etamax = 0;

%spline order
PAR.c = 4; 
PAR.L1 = 50; %# of steps for body along length
PAR.L2 = 20; %# of steps for head along length
PAR.L3 = 30; %# of steps for wing around the boundary
PAR.T1 = 25; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;

%Load initial condition
load([PAR.solutionpath 'fly_' PAR.stub '000' '/' 'ManualFit_' PAR.stub '000']);

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;

% ------------------------------------
% Load images
for cam=1:PAR.numcam
    input_filename = sprintf('%s/%s/cam%03d/%s%06d.bmp', ...
        PAR.videopath,PAR.filetag,cam,PAR.filetag,frame);
    Input(:,:,cam) = imread(input_filename);
end

PAR.imgres = size(Input(:,:,1));
PAR.imgres_crop = PAR.imgres;

% ------------------------------------
% Load the solution data
load([PAR.solutionpath 'fly_' PAR.stub '000' '/fly' num2str(frame) '.mat']);
load([PAR.solutionpath 'Features_' PAR.stub '000' '/Features' num2str(frame) ...
    '.mat']);

% load([PAR.solutionpath 'fly_' PAR.stub '_BodyPredict/fly' num2str(frame) '.mat']);
% load([PAR.solutionpath 'Features_' PAR.stub '_BodyPredict/Features' num2str(frame) ...
%     '.mat']);

%----------------------------
% Evaluate Model at solution
%----------------------------
% Manual fit solution at first frame
% sol = reshape(ManualFit.soln,[],1);

% Predicted solution at current frame
%sol1 = InternalVariablesDS.xh_';

% Actual solution
sol = xh;

[x,y,z] = flymodQ(sol,PAR.params,PAR);

Flypts = cell(3,1);
R = [];
for n = 1:length(x)
    Flypts{n} = [reshape(x{n},[],1) reshape(y{n},[],1) reshape(z{n},[],1)];
    C(:,n) = PAR.cam(n).C;
    R = [R ; PAR.cam(n).R];
end
pts = Flypts;
Flypts = cell2mat(Flypts)';
Flypts = [Flypts ; zeros(1,size(Flypts,2))];

%drawscene(Flypts,C,R,fignum,'none','Fly Experimental Setup');
% %view(0,0);

%Calcualte the projection Rays
%Features = feat_detect(frame,[],PAR);

figure(fignum);
hold on;

%Now, draw them
linelen = [140 120 130];
crs = {'c','y','m'};
for cam = 1%:PAR.numcam
    Dirs = Features(frame).DataRays{cam}(:,1:3);

    %Get the points that are corresponding to the 2D boundary
    pts3D = [];
    %iterate over each part
    for i = 1:length(pts)
        index = Features(frame).IdxTo3DPts{cam}{i};
        temp = pts{i}(index,:);
        pts3D = [pts3D ; temp];
    end

    % Get rid of the points that are occluded
    occidx = Features(frame).occluded_idx{cam};
    pts3D(occidx,:) = [];

    npts = 2;
    for j = 1:size(Dirs,1)
        %Calculate the scalar distance lambda from the camera center to the
        %3D point.  x = C + lambda*Dirs(j,:);
        dd = (pts3D(j,:)-PAR.cam(cam).C') * Dirs(j,:)';
        
        linepts = repmat(PAR.cam(cam).C,1,npts) + ...
            repmat(linspace(0,dd,npts),3,1).*repmat(Dirs(j,:)',1,npts);
        plot3(linepts(1,:),linepts(2,:),linepts(3,:),'-','color',crs{cam},'linewidth',.5)
    end

%     npts = 5;
%     for j = 1:size(Dirs,1)
%         linepts = repmat(PAR.cam(cam).C,1,npts) + ...
%             repmat(linspace(0,linelen(cam),npts),3,1).*repmat(Dirs(j,:)',1,npts);
%         plot3(linepts(1,:),linepts(2,:),linepts(3,:),'-','color',crs{cam},'linewidth',.5)
%     end
end

% Now Draw Fly surface
colors = {'g','g','g'};
for i = 1:length(x)
    if i > 1
        color = [.8 .8 .8];
    else
        color = [.314 .314 .314];
    end
    surf(x{i},y{i},z{i},'facecolor',color,'facelighting','phong',...
        'specularstrength',0.1);
end

% set(gca,'zdir','reverse','ydir','reverse',...
%     'xlim',[10 15],'ylim',[3 9],'zlim',[13 17],...
%     'xticklabel',[],'yticklabel',[],'zticklabel',[]);
set(gca,'zdir','reverse','ydir','reverse',...
    'xticklabel',[],'yticklabel',[],'zticklabel',[]);
%axis equal
view(-77,7.26)
grid off;
%box;
camlight;
axis tight

% %closeup View of point correspondence
% for cam = 1:PAR.numcam
%     
%     figure(fignum*10+cam);
%     for i = 1:length(x)
%         surf(x{i},y{i},z{i},'facecolor',colors{i});
%         hold on;
%     end
%     
%     axis equal
%     pts3D = [];
%     %iterate over each part
%     for i = 1:length(pts)
%         index = Features(frame).IdxTo3DPts{cam,1}{i};
%         temp = pts{i}(index,:);
%         pts3D = [pts3D ; temp];
%     end
%     
%     % Get rid of the points that are occluded
%     occidx = Features(frame).occluded_idx{cam,1};
%     pts3D(occidx,:) = [];
%     
%     %line directions: First 3 components of ray's Plüker
%     %coordinates
%     Direction = Features(frame).DataRays{cam}(:,1:3);
%     
%     for ii = 1:size(pts3D,1)
%         plot3(pts3D(ii,1),pts3D(ii,2),pts3D(ii,3),'yo');
%         npts = 500;
%         
%         linepts = repmat(PAR.cam(cam).C,1,npts) + repmat(linspace(0,200,npts),3,1).*repmat(Direction(ii,:)',1,npts);
%         
%         plot3(linepts(1,:),linepts(2,:),linepts(3,:),'c-')
%         
%         cp = kdtree(linepts',pts3D(ii,:));
%         dist{cam}(ii,:) = sqrt(sum( diff([pts3D(ii,:) ; cp],1,1).^2,2));
%         
%         plot3([pts3D(ii,1) cp(1)],[pts3D(ii,2) cp(2)],[pts3D(ii,3) cp(3)],'r-')
%     end
% end
