%MakeSynVideo.m
%
% this function makes a set of sequentially numbered bitmaps that are
% rendered movies of a solution that has already been tracked (i.e. a
% synthetic video).  
%
% It saves then in the specified folder

clear all
close all

%Load the previously stored ManualFit Parameters
[FileName,PathName] = uigetfile({'*.mat'},'Select "ManualFit" data file for the video sequence');
load([PathName FileName]);

PAR = ManualFit.ImageData;

fignum = 1;

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .8;

% Model color
color = {'r','r'};
lw = .5;


%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
frames = [25 135];
%frames = ManualFit.ImageData.frames4analysis; % number of frames in movie

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


% Make directories to save the synthetic video if they don't exist
curdir = pwd;
cd(PAR.imagepath);
cd ..
ImageFolder = pwd;
if exist([PAR.stub '_syn'],'dir') ~= 7
    mkdir([PAR.stub '_syn']);
    %make camera directories too
    for i = 1:PAR.numcam
        mkdir([PAR.stub '_syn'],['cam00' num2str(i)]);
    end
end
cd(curdir);


soln1 = zeros(PAR.numfly*PAR.statedim,length(movidx));

if interptag == 1
    for i=1:length(movidx)
        load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(movidx(i)) ...
            '.mat']);
        soln1(:,i) = xh(1:PAR.numfly*PAR.statedim);
        clear xh InternalVariablesDS
    end
end

%Now, perform smoothing of the tracked data
for k = 1:size(soln1,1)
    if k <=7
        [sp,soln1(k,:)] = spaps(frames(1):frames(2),soln1(k,:),.05);
    else
        [sp,soln1(k,:)] = spaps(frames(1):frames(2),soln1(k,:),.5);
    end
end

%renormalize
iidx = 4:7;
Mag = sqrt(sum(soln1(iidx,:).^2,1));
soln1(iidx,:) = soln1(iidx,:)./repmat(Mag,4,1);

iidx = 8:11;
Mag = sqrt(sum(soln1(iidx,:).^2,1));
soln1(iidx,:) = soln1(iidx,:)./repmat(Mag,4,1);

iidx = 12:15;
Mag = sqrt(sum(soln1(iidx,:).^2,1));
soln1(iidx,:) = soln1(iidx,:)./repmat(Mag,4,1);



%%
for i= 1:size(soln1,2)
    kk = frames(1)+(i-1);

    
    clear flymod
    [x,y,z] = flymodQ(soln1(:,i),PAR.params,PAR);
    [x1,y1,z1] = flymodQ_WingThick(soln1(:,i),PAR.params,PAR);

    %----------------------
    % Render the model to an image
    %----------------------
    pts = [];
    for j = 1:length(x)
        PAR.modsample(j) = size(x{j},1);
        PAR.modsample1(j) = size(x1{j},1);
        pts{j} = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
        pts1{j} = [reshape(x1{j},[],1) reshape(y1{j},[],1) reshape(z1{j},[],1)];
    end
    
    for m = 1:PAR.numcam %Iterate over each camera
        %Return camera parameters
        R = PAR.cam(m).R;
        K = PAR.cam(m).K;
        camC = PAR.cam(m).C;

        %First, render the 3D model points to the current camera view

        % Return only the image of the fly
        IMout = renderflySYN(pts,PAR,m,pts1);
        IMout = uint8(IMout);
        
        %keyboard

        % ------------------------------------
        % write the image to disk
        digits = length(num2str(PAR.numframes));
        image2write = sprintf(['%s%0' num2str(digits) 'd%s'],[PAR.stub '_syn'],kk,PAR.image_filter(2:end));
        path2write = sprintf('%s/%s/cam%03d/',ImageFolder,[PAR.stub '_syn'],m);
        imwrite(IMout,[path2write image2write],PAR.image_filter(3:end));

    end
end

