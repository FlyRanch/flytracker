% Paste_image.m
%
% This function plots the tracking results for a particular frame of a
% movie sequence defined by the user
%
% The result is created by plotting the image in one axis (ax1) and 
% plotting the model result and feature points in an identical axis (ax2)
% on top of the image.
%
% parameters at the beginning of the file need to be changed accordingly
%================================
%
% If a file "ManualFit" already exists in memory, this function will use it
% to plot the data.
%
% If not, it will query the user to select a sequence.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all; 
clc;
if exist('ManualFit','var') == 1
    PAR = ManualFit.ImageData;
else
    [FileName,PathName] = uigetfile({'*.mat'},'Select "ManualFit" data file for the video sequence');
    load([PathName FileName]);
       
    % Check that the paths stored in 'ManualFit' match the location that you
    % just selected the file from
    % If 75% of the paths match , then we'll assume everything's okay.
    endd = round(.75*length(ManualFit.ImageData.solutionpath));
    
    if strcmp(PathName(1:endd),ManualFit.ImageData.solutionpath(1:endd)) == 0
        %if different, run Loadvideo routine
        warning('The directories stored in ManualFit structure do not match its current location.\n You are prompted to relocate the directories',[]);
        PAR = LoadVideo;
        ManualFit.ImageData = PAR;
        save([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub],'ManualFit');
    else
        PAR = ManualFit.ImageData;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure number to plot in and frame to plot
fignum = 1;
frame = 420;
for k = 1:3
   close(10*fignum+k)
end

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = 1;

% Model color
color = {'g','g'};
lw = 1;

plotpredict = 1;
plotcur = 1;
OccTag = 0;    %Set to 1 to see the artificial occlusion

cptsTag = 1;   %Set to 1 to see the corresponding points
KF_tag = 0;    %Set to 1 to see the 1D measurement lines searching for edge points
nrmlTag = 0;   %Set to 1 to see the normal vectors
renIMtag = 0; 

PAR.pixpermm = 1;
PAR.numfly = 1;
%Number of parameters of the model (i.e. 8 control points)
PAR.mdlpar = 15*ones(1,PAR.numfly);
PAR.modelfun_H = @modcurvesplineP;
PAR.etamax = 0;

%spline order
PAR.c = 4; 
PAR.L1 = 15; %# of steps for body along length
PAR.L2 = 6; %# of steps for head along length
PAR.L3 = 20;%30; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing
PAR.T3 = 3; %# of steps along thickness of wing

% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;

% Create an artificial occlusion in each camera view.  Each cell represents
% the vertices of a polygon in the image for the particular camera
PAR.OccludeShape = {RectVertices(500,25,256,230),...
    RectVertices(500,20,256,240),[]}; 

%Load initial condition
load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;

% ------------------------------------
% Load images
clear Input
for cam=1:PAR.numcam
    %Number of digit places for the total number of frames
    digits = length(num2str(PAR.numframes));

    input_filename = sprintf(['%s/cam%03d/%s%0' num2str(digits) 'd%s'], ...
        PAR.imagepath, cam, PAR.stub, frame,PAR.image_filter(2:end));

    Input(:,:,cam) = imread(input_filename);
end

PAR.imgres = size(Input(:,:,1));
PAR.imgres_crop = PAR.imgres;

% ------------------------------------
% % Load the solution data
load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(frame) '.mat']);
load([PAR.solutionpath 'Features_' PAR.stub '/Features' num2str(frame) ...
    '.mat']);

% load([PAR.solutionpath 'fly_' PAR.stub '_FD_BodyScale/fly' num2str(frame) '.mat']);
% load([PAR.solutionpath 'Features_' PAR.stub '_FD_BodyScale/Features' num2str(frame) ...
%     '.mat']);

%----------------------------
% Evaluate Model at solution
%----------------------------
% Manual fit solution at first frame
% sol = reshape(ManualFit.solnQ,[],1);

% load testSigpts
% sol = testSigpts(:,10);

% Actual solution
% load exp101_ManualTrack;
% sol = State(:,frame);
% sol = InternalVariablesDS.xh_;
sol = xh;

if plotpredict == 1
    load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(frame+1) '.mat']);
    sol1 = InternalVariablesDS.xh_;
    
%     bodyscrew = sol1(1:6);
%     Gbody = screw2homo(bodyscrew);
% 
%     Tbody = Gbody(1:3,4);
%     Rbody = Gbody(1:3,1:3);
%     q_body = quat2matNEW(Rbody);
%     qL = [0.3028
%         0.3734
%         0.1815
%         0.8579];
%     qR = [-0.2691
%         0.5826
%         -0.0903
%         0.7616];
%     sol1 = [Tbody;q_body;qL;qR];
    
    load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(frame) '.mat']);
else
    % Predicted solution at current frame
    sol1 = InternalVariablesDS.xh_;
end
    
clear flymodQ
%clear flymod_BodyOnly
[x,y,z] = flymodQ(sol,PAR.params,PAR);
% [x1,y1,z1] = flymodQ(sol1,PAR.params,PAR);
[x1,y1,z1] = flymodQ(sol1,PAR.params,PAR);

%----------------------
% Plot Movie Frame or Rendered model image
%----------------------
pts = [];
for i = 1:length(x)
    PAR.modsample(i) = size(x{i},1);
    PAR.modsample1(i) = size(x1{i},1);
    
    pts{i} = [reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
    if i == 1
        % get dorsal edge of the body
        dorsalpts = [x{i}(:,10),y{i}(:,10),z{i}(:,10)];
        %dorsalpts = [x{i}(:,6),y{i}(:,6),z{i}(:,6)];
%         dorsalpts = [x{i}(end,4),y{i}(end,4),z{i}(end,4)
%             x{i}(8,4),y{i}(8,4),z{i}(8,4)
%             x{i}(8,7),y{i}(8,7),z{i}(8,7)];
    end
    pts1{i} = [reshape(x1{i},[],1) reshape(y1{i},[],1) reshape(z1{i},[],1)];
end

% Initialize projected planar points
u = cell(PAR.numcam,length(x));
v = cell(PAR.numcam,length(x));
uT = cell(PAR.numcam,length(x));
vT = cell(PAR.numcam,length(x));

u1 = cell(PAR.numcam,length(x));
v1 = cell(PAR.numcam,length(x));
uT1 = cell(PAR.numcam,length(x));
vT1 = cell(PAR.numcam,length(x));


for i = 1:PAR.numcam
        
    DLTparams = PAR.DLT(:,i);
    
    shiftax = [-1 -1 -1];
    xax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 11 5 14]+ repmat(shiftax,2,1));
    yax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 6 14]+ repmat(shiftax,2,1));
    zax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 5 15]+ repmat(shiftax,2,1));
    
    figure(fignum*10+i);
    
    imagesc(Input(:,:,i));
    %imagesc(imadjust(Input(:,:,i)));
    colormap gray

    ax1 = gca;
    %set(ax1,'visible','off');
    set(ax1,'units','pixels','visible','off','position',[0 0 PAR.imgres(2) PAR.imgres(1)].*plotscale);
    axis tight
    
    aspectratio = get(ax1,'PlotBoxAspectRatio');
    impos = get(ax1,'position');
    ax2 = axes('position',[0 0 .5 .5],'PlotBoxAspectRatio',aspectratio);
    
    subax(2*i-1) = ax1;
    subax(2*i) = ax2;
    
    for j = 1:length(pts)
        
        dorsal2Dpts = dlt_3D_to_2D(DLTparams,dorsalpts);
        uv = dlt_3D_to_2D(DLTparams,pts{j});
        uv1 = dlt_3D_to_2D(DLTparams,pts1{j});
        
        
        u{i,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2));
        v{i,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2));
        u1{i,j} = reshape(uv1(:,1),size(x1{j},1),size(x1{j},2));
        v1{i,j} = reshape(uv1(:,2),size(x1{j},1),size(x1{j},2));
        
        uT{i,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2))';
        vT{i,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2))';
        uT1{i,j} = reshape(uv1(:,1),size(x1{j},1),size(x1{j},2))';
        vT1{i,j} = reshape(uv1(:,2),size(x1{j},1),size(x1{j},2))';
        
        if j == 1
            %Plot dorsal edge
            plot(ax2,dorsal2Dpts(:,1),dorsal2Dpts(:,2),'r-','linewidth',5);
            
            hold on;
            
            %Plot Head and tail
            %plot(ax2,uv([1 end],1),uv([1 end],2),'rs');
        end
        
        
        if plotcur == 1 
            patch(u{i,j},v{i,j},[1 1 1],'facecolor','none','edgecolor','r','linewidth',lw);
            patch(uT{i,j},vT{i,j},[1 1 1],'facecolor','none','edgecolor','r','linewidth',lw);
            %patch(u{i,j},v{i,j},[1 1 1],'facecolor',[.8 .8 .8],'edgecolor','w');
            %patch(uT{i,j},vT{i,j},[1 1 1],'facecolor',[.8 .8 .8],'edgecolor','w');
        end
        if plotpredict == 1 
            patch(u1{i,j},v1{i,j},[1 1 1],'facecolor','none','edgecolor','y');
            patch(uT1{i,j},vT1{i,j},[1 1 1],'facecolor','none','edgecolor','y');
        end
        
        
    end
    
    if (OccTag && ~isempty(PAR.OccludeShape{i}) )
        Vert = PAR.OccludeShape{i};
        patch(Vert(:,1),Vert(:,2),[1 1 1],'facecolor','none','edgecolor','r');
    end
    
    plot(xax(:,1),xax(:,2),'r.-',yax(:,1),yax(:,2),'g.-',zax(:,1),zax(:,2),'b.-')
    %label the axes
    text(xax(2,1),xax(2,2),'X','fontsize',12)
    text(yax(2,1),yax(2,2),'Y','fontsize',12)
    text(zax(2,1),zax(2,2),'Z','fontsize',12)

    %--- Plot Correspondence between points
    if cptsTag == 1;
%         datapts = Features(frame).DataptsFullIC{i};
%         mdlpts = Features(frame).ModelptsIC{i};
%         
        datapts = Features(frame).DataptsFull{i};
        mdlpts =  Features(frame).Modelpts{i};
        sz = 5;
        ww = 1;
        for j = 1:size(datapts,1)
            if datapts(j,:) == [0 0]
                %This point is occluded. Plot as yellow star.
                plot(mdlpts(j,1),mdlpts(j,2),'y*')
            else
                %plot([mdlpts(j,1) datapts(j,1)],[mdlpts(j,2) datapts(j,2)], ...
                %    ['co:'],'markerfacecolor','k');
                         plot([mdlpts(j,1) datapts(j,1)],[mdlpts(j,2) datapts(j,2)],'c--','linewidth',ww);
                          plot(mdlpts(j,1),mdlpts(j,2),'co','markerfacecolor','k','markersize',sz);
                          plot(datapts(j,1),datapts(j,2),'ko','markerfacecolor','c','markersize',sz);
            end
        end

        %--------------
        % Plot all the boundary points and high curvature points
        %--------------
        %     plot(Features(frame).All_Bndy_Points(:,1), ...
        %          Features(frame).All_Bndy_Points(:,2),'b.');
        %
        %     plot(Features(frame).Candpts(:,1), ...
        %          Features(frame).Candpts(:,2),'rs','markersize',10);

    end
    
    if nrmlTag == 1
        
        %--------------------------------------------
        % Plot all the outward normal vectors from the projected model
        % contour.  These directions are used for edge feature detection.
        %---------------------------------------------
        quiver(Features(frame).Modelpts{i}(:,1), ...
            Features(frame).Modelpts{i}(:,2),Features(frame).Nrml{i}(:,1),...
            Features(frame).Nrml{i}(:,2),1.0,'y')
    end
    
    %=====================================================================
    if  KF_tag == 1
        %Grab the Normal Vectors
        Nrml = Features(frame).Nrml{i};
        PredBdypts = Features(frame).Modelpts{i};
        
        for n = 1:size(PredBdypts,1)
            r = n;
            NN = Nrml(r,:);
            RR = [15 15];% sqrt(NRy(n));
            lambda = linspace(-RR(1),RR(2),30);
            X = repmat(PredBdypts(n,:),length(lambda),1) + repmat(lambda,2,1)'.*repmat(NN,length(lambda),1);
            hold on
            plot(X(:,1),X(:,2),'r-','linewidth',1);


            %            if ~isempty(Features(frame).Edgepts{r,k})
            %                plot(Features(frame).Edgepts{r,k}(:,1),Features(frame).Edgepts{r,k}(:,2),'o','color',color{k})
            %            end
        end
    end

    set(ax2,'units','pixels','fontsize',12,'position',impos,'color','none','xlim',...
        [0.5 PAR.imgres(2)+0.5],'ylim',[-0.5 PAR.imgres(1)-0.5],'visible','off',...
        'xdir','normal','ydir','reverse');
    
    %set(fignum,'position',[0 500 PAR.imgres(2)*plotscale PAR.imgres(1)*plotscale])
    set(fignum*10+i,'units','pixels','position',[20+(i-1)*PAR.imgres(2) 50 PAR.imgres(2) PAR.imgres(1)].*plotscale);
    %title(['Frame ' num2str(frame) ', IC: p0 = ' num2str([sol(1:8)]) ])
    text(PAR.imgres(2)/2 - 50,20,['\color{white}Frame ' num2str(frame)], ...
        'fontsize',20);
end

if renIMtag == 1
    for i = 1:PAR.numcam

        %Render Model
        [Y,idx,IMout,Nrml] = renderflyMOD2(pts,PAR,i,pts1);

        figure;
        imagesc(IMout); colormap gray;
        hold on;
        for j = 1:length(Y)
            axes(subax(2*i));
            plot(Y{j}(:,1),Y{j}(:,2),'b.');
            
            quiver(Y{j}(:,1),Y{j}(:,2),Nrml{j}(:,1),Nrml{j}(:,2),'r')
        end
    end
end