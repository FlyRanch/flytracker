close all;clear all;clc;


global PAR
addpath('mex/');
addpath('core/');

startpath = '/Users/cisco/FlyVideo/solutions/';
getImageData = true;

PAR = LoadVideo;

ManualFit.ImageData = PAR;

% Make two directories to save the estimated state and the features into
% if they don't already exist
if exist([PAR.solutionpath ['fly_' PAR.stub]],'dir') ~= 7
    mkdir(PAR.solutionpath,['fly_' PAR.stub]);
    mkdir(PAR.solutionpath,['Features_' PAR.stub]);
end

save([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub],'ManualFit');

%exp98: 325:536
startframe = 326;
endframe = 536;%PAR.numframes;
numframes = 210;
PAR.frames2skip = [];
PAR.framesample = 1;

% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;
PAR.imgres = [512 512];
%PAR.CalibFile = 'kine/Calibrations for Kine/cal_for_20031201_done_20060203.mat';
%PAR.CalibFile = 'kine/Calibrations for Kine/cal_for_20031202_done_20060127.mat';
%PAR.CalibFile = 'kine/Calibrations for Kine/cal_of_20031208_on_20040129';
PAR.CalibFile = 'kine/Calibrations for Kine/cal_for_20031217_done_20040205';

%exp098: 
PAR.ICframe = 326;
PAR.FileFromKine = 'cisco098';

getIC = 1;
getBodyShape = 1;

%% ====================================================
% START ADVANCED PARAMETERS
%=====================================================

PAR.flynum = [1 ];  %[1 2];
PAR.numfly = length(PAR.flynum);
PAR.OccludeShape = {[],[],[]};
%spline order
PAR.c = 4; 
PAR.L1 = 15; %# of steps for body along length
PAR.L2 = 6; %# of steps for head along length
PAR.L3 = 30; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

PAR.etamax = 0;
PAR.paramdim = 1;
PAR.statedim = 15*ones(1,PAR.numfly); 
PAR.pNoisedim = PAR.statedim;

PAR.NumPrevSol = 5;

PAR.streams = 2;

% Dimension of the parameters that define the model shape
% for fly its [R L]
% R = 20 control points that define the B-spline function for the
% radius function R(u); size(R) = [1 20];
% L = [Length of tail region , Length of body/head region]; size(L) = [1 2];
% PAR.paramdim = 22;
% PAR.paramdim = 10; (for fish03)

% END ADVANCED PARAMETERS
frame = 326:520;
%=======================================================
%% Get Initial Conditions
for ii = 1
    ManualFit.ImageData = PAR;

    %Calculate the initial condition
    ManualFit = auto_init(ManualFit,frame(ii),'BG');
    ManualFit = auto_init(ManualFit,frame(ii),'IC');
    save([PAR.solutionpath 'fly_' PAR.stub '/ManualFit_' PAR.stub],'ManualFit');
    
    clear ManualFit
    display('%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause
end