%plot_syn_error_metrics.m
%
%plot_syn_error_metrics plots the state variables of a specified sequence over a
%specified number of frames.
%It will plot the estimated values from traking the synthetic images and the 
%known values that the synthetic images were created from.
%
clear all
close all

PAR.videopath = '../FlyVideo/';

%The first one is the tracked states, and the second is the known states
%used to create the synthetic images.
PAR.filetag = {'exp098000_syn','exp098000'};

PAR.solutionpath = [PAR.videopath '/solutions/'];
frgroup = {[325 521],[325 521]};

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

PAR.stub = [PAR.filetag{1}];
fignum = 1;

% --------------------------------
% Load the ManualFit Parameters
load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);

%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
frames = frgroup{1}; % number of frames in movie

samplerate = 1;
movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);

solidx = [1 frames(2)];


% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;


soln1 = zeros(PAR.numfly*PAR.statedim,length(movidx));
soln2 = zeros(PAR.numfly*PAR.statedim,length(movidx));


for i=1:length(movidx)
    load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(movidx(i)) ...
        '.mat']);
    soln1(:,i) = xh(1:PAR.numfly*PAR.statedim);

    load([PAR.solutionpath 'fly_' PAR.filetag{2} '/fly' num2str(movidx(i)) ...
        '.mat']);
    soln2(:,i) = xh(1:PAR.numfly*PAR.statedim);
    clear xh InternalVariablesDS
end

%Now, perform smoothing of the tracked data
ii = 5;
for k = 1:size(soln2,1)
    if k <=7
        [sp,soln2(k,:)] = spaps(frames(1):frames(2),soln2(k,:),.05);
    else
        [sp,soln2(k,:)] = spaps(frames(1):frames(2),soln2(k,:),.5);
    end
end

%renormalize
iidx = 4:7;
Mag = sqrt(sum(soln2(iidx,:).^2,1));
soln2(iidx,:) = soln2(iidx,:)./repmat(Mag,4,1);

iidx = 8:11;
Mag = sqrt(sum(soln2(iidx,:).^2,1));
soln2(iidx,:) = soln2(iidx,:)./repmat(Mag,4,1);

iidx = 12:15;
Mag = sqrt(sum(soln2(iidx,:).^2,1));
soln2(iidx,:) = soln2(iidx,:)./repmat(Mag,4,1);

%% ===============================================================

% Convert the automatically captured body data into Euler components
BodyTrans_auto = zeros(3,length(movidx));
BodyAng_auto = zeros(3,length(movidx));

BodyTrans = zeros(3,length(movidx));
BodyAng = zeros(3,length(movidx));

for k = 1:length(movidx)
    Rbod = quat2matNEW(soln1(4:7,k));
    BodyAng_auto(:,k) = Rot2Euler(Rbod');%Rot2Euler(G(1:3,1:3)');
    BodyTrans_auto(:,k) = soln1(1:3,k);%G(1:3,4);

    Rbod = quat2matNEW(soln2(4:7,k));
    BodyAng(:,k) = Rot2Euler(Rbod');%Rot2Euler(G(1:3,1:3)');
    BodyTrans(:,k) = soln2(1:3,k);%G(1:3,4);
end

%Common frames between manual and automatically tracked data
ff = 1:length(movidx); %frames(1):frames(2);
BodyAng = BodyAng.*(180/pi);
BodyAng_auto = BodyAng_auto.*(180/pi);
tt = (0:length(ff)-1)*PAR.dt*1000;

figure;
subplot(1,2,1);
mkr = {'--','-','o','*'};
hold on;


%plot prior database
hh = plot(tt,BodyAng(1,ff),['r' mkr{1}],tt,BodyAng(2,ff),['g' mkr{1}],tt,BodyAng(3,ff),['b' mkr{1}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',1)
end


%plot tracked data
hh = plot(tt,BodyAng_auto(1,ff),['r' mkr{2}],tt,BodyAng_auto(2,ff),['g' mkr{2}],tt,BodyAng_auto(3,ff),['b' mkr{2}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',2)
end

%set(gca,'ytick',-pi:pi/4:pi,'xlim',[ff(1) ff(end)]);
set(gca,'xlim',[tt(1) tt(end)]);
%legend(hh,'roll','pitch','yaw');
title(PAR.stub);
ylabel('degrees');
xlabel('ms');
grid

subplot(1,2,2);
hold on;

%plot prior database
hh = plot(tt,BodyTrans(1,ff),['r' mkr{1}],tt,BodyTrans(2,ff),['g' mkr{1}],tt,BodyTrans(3,ff),['b' mkr{1}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',1)
end


%plot tracked data
hh = plot(tt,BodyTrans_auto(1,ff),['r' mkr{2}],tt,BodyTrans_auto(2,ff),['g' mkr{2}],tt,BodyTrans_auto(3,ff),['b' mkr{2}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',2)
end

set(gca,'xlim',[tt(1) tt(end)])
%legend(hh,'X','Y','Z');
title('Body Translation');
ylabel('mm');
xlabel('ms');
grid
%keyboard


%Normalize by body length
BodyTrans = BodyTrans ./ (ManualFit.params.bodyscale*(...
    ManualFit.params.bodylen + ManualFit.params.headlen));

BodyTrans_auto = BodyTrans_auto ./ (ManualFit.params.bodyscale*(...
    ManualFit.params.bodylen + ManualFit.params.headlen));

%% Error in body transformation
AngErr = (BodyAng_auto(:,ff) - BodyAng(:,ff));


%% Wing Coordinates

LwingAng_auto = zeros(3,length(movidx));
RwingAng_auto = zeros(3,length(movidx));

LwingAng = zeros(3,length(movidx));
RwingAng = zeros(3,length(movidx));

for k = 1:length(movidx)
    wing = quat2matNEW(soln1(8:11,k));
    LwingAng_auto(:,k) = Rot2Euler(wing).*(180/pi);
    
    wing = quat2matNEW(soln1(12:15,k));
    RwingAng_auto(:,k) = Rot2Euler(wing).*(180/pi);
    
    wing = quat2matNEW(soln2(8:11,k));
    LwingAng(:,k) = Rot2Euler(wing).*(180/pi);
    
    wing = quat2matNEW(soln2(12:15,k));
    RwingAng(:,k) = Rot2Euler(wing).*(180/pi);
end


figure; 
% --------------------------------
% Left wing
subplot(2,1,1)
%plot prior database
hh = plot(tt,LwingAng(1,ff),['r' mkr{1}],tt,LwingAng(2,ff),['g' mkr{1}],tt,LwingAng(3,ff),['b' mkr{1}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',1)
end


%plot tracked data
hh = plot(tt,LwingAng_auto(1,ff),['r' mkr{2}],tt,LwingAng_auto(2,ff),['g' mkr{2}],tt,LwingAng_auto(3,ff),['b' mkr{2}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',2)
end

%set(gca,'ytick',-pi:pi/4:pi,'xlim',[ff(1) ff(end)]);
set(gca,'xlim',[tt(1) tt(end)]);
%legend(hh,'roll','pitch','yaw');
title(PAR.stub);
ylabel('degrees');
xlabel('ms');
grid

% --------------------------------
% Right wing
subplot(2,1,2);
hh = plot(tt,RwingAng(1,ff),['r' mkr{1}],tt,RwingAng(2,ff),['g' mkr{1}],tt,RwingAng(3,ff),['b' mkr{1}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',1)
end


%plot tracked data
hh = plot(tt,RwingAng_auto(1,ff),['r' mkr{2}],tt,RwingAng_auto(2,ff),['g' mkr{2}],tt,RwingAng_auto(3,ff),['b' mkr{2}]);

for i = 1:length(hh)
    set(hh(i),'linewidth',2)
end

%set(gca,'ytick',-pi:pi/4:pi,'xlim',[ff(1) ff(end)]);
set(gca,'xlim',[tt(1) tt(end)]);
%legend(hh,'roll','pitch','yaw');
title(PAR.stub);
ylabel('degrees');
xlabel('ms');
grid

%% --------------------------------
% %Quaternion values
% % Left wing
% subplot(2,1,1);
% %plot prior database
% hh = plot(ff,soln2(8,ff),['r' mkr{1}],ff,soln2(9,ff),['g' mkr{1}],ff,soln2(10,ff),['b' mkr{1}],ff,soln2(11,ff),['k' mkr{1}]);
% 
% hold on;
% %plot tracked data
% hh = plot(ff,soln1(8,ff),['r' mkr{2}],ff,soln1(9,ff),['g' mkr{2}],...
%     ff,soln1(10,ff),['b' mkr{2}],ff,soln1(11,ff),['k' mkr{2}]);
% 
% %set(gca,'ytick',-pi:pi/4:pi);
% %legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
% legend(hh,'q_1','q_2','q_3','q_4');
% title('Left wing');
% grid
% 
% % --------------------------------
% % Right wing
% subplot(2,1,2);
% %plot Known states
% hh = plot(ff,soln2(12,ff),['r' mkr{1}],ff,soln2(13,ff),['g' mkr{1}],ff,soln2(14,ff),['b' mkr{1}],ff,soln2(15,ff),['k' mkr{1}]);
% 
% hold on;
% %plot tracked data
% hh = plot(ff,soln1(12,ff),['r' mkr{2}],ff,soln1(13,ff),['g' mkr{2}],...
%     ff,soln1(14,ff),['b' mkr{2}],ff,soln1(15,ff),['k' mkr{2}]);
% 
% % set(gca,'ytick',-pi:pi/4:pi);
% % legend(hh,'R_x - roll','R_y - pitch','R_z - yaw');
% legend(hh,'q_1','q_2','q_3','q_4');
% title('Right wing');
% grid