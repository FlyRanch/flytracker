%show_sol_3Dpath.m
%
%show_sol_auto makes a movie of the tracking result for a particular video
%sequence.  It works similar to 'paste_image_auto' but iterates over each
%frame
%
%After running this program the user should call 'movie2avi(M2,...)' to
%save the movie to whatever filename they choose.

clear all;close all;clc

PAR.videopath = '../FlyVideo/';
PAR.filetag = 'exp095000';

% PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.solutionpath = ['/solutions/'];
PAR.stub = [PAR.filetag];
fignum = 1;

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .8;


% --------------------------------
% Load the ManualFit Parameters
% load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);
load([cd PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub '.mat']);

%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
frames = [1 120]; % number of frames in movie

samplerate = 1;
movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);



cptsTag = 0;     %Set to 1 to see the corresponding points
nrmlTag = 0;     %Set to 1 to see the normal vectors
KF_tag = 0;      %Set to 1 to see the Predicted Model and Gating Ellipses
interptag = 1;
%--------------------------------
% Define the Tracking Parameters
PAR.pixpermm = 1;
PAR.numfly = 1;
%Number of parameters of the model (i.e. 8 control points)
PAR.mdlpar = 15*ones(1,PAR.numfly);
PAR.statedim = PAR.mdlpar;
PAR.modelfun_H = @modcurvesplineP;
PAR.etamax = 0;

%spline order
PAR.c = 4;
PAR.L1 = 15; %# of steps for body along length
PAR.L2 = 6; %# of steps for head along length
PAR.L3 = 25; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;
solidx = [1 length(movidx)];

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;


soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
SOLN = zeros(length(frames(1):frames(2)),PAR.numfly*PAR.statedim);
if interptag == 1
    for i=1:length(movidx)
        load([cd PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(movidx(i)+430-1) ...
            '.mat']);
        soln1(i,:) = xh(1:PAR.numfly*PAR.statedim)';
        clear xh InternalVariablesDS
    end
end

%Now, perform interpolation to calculate the state at in between frames
for k = 1:size(soln1,2)
    SOLN(:,k) = interp1(movidx,soln1(:,k),frames(1):frames(2),'spline')';
end

%% Calculate the region of fly motion
buffer = 1;
minx = min(SOLN(:,1)) - buffer;
maxx = max(SOLN(:,1)) + buffer;
miny = min(SOLN(:,2)) - buffer;
maxy = max(SOLN(:,2)) + buffer;
minz = min(SOLN(:,3)) - buffer;
maxz = max(SOLN(:,3)) + buffer;

%--------------------------------------------
%calculate theta-dot
theta_vec = sqrt(sum(SOLN(:,4:6).^2,2));
theta_dot = gradient(theta_vec,PAR.dt);
dt = PAR.dt;
tsteps = 1;
dummy = [];
torder = 1;

M1dum = [];
M2dum = [];

T = 0;%linspace(0,12,60);

M2 = moviein(size(SOLN));
for i= 1:size(SOLN,1)
    kk = frames(1)+(i-1);

    clear flymodQ
    [x,y,z] = flymodQ(SOLN(i,:),PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

       
    
    
    figure(fignum);
    clf;
    
    for k = 1:length(x);
        surf(x{k},y{k},z{k},'facecolor','r','edgecolor','k','facelighting','phong');
        hold on;
    end

    set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz],...
        'zdir','reverse','ydir','reverse');
    %axis vis3d
    
    figure(fignum);
    set(fignum,'color','w');
    %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])
    

    M2(i) = getframe(fignum);
%     M2(kk) = getframe(fignum);
end

