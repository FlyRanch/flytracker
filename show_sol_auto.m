%show_sol_auto.m
%
%show_sol_auto makes a movie of the tracking result for a particular video
%sequence.  It works similar to 'paste_image_auto' but iterates over each
%frame
%
%After running this program the user should call 'movie2avi(M2,...)' to
%save the movie to whatever filename they choose.

clear all;close all;clc

PAR.videopath = '../FlyVideo/';
PAR.filetag = 'exp101000';

PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.stub = [PAR.filetag];
fignum = 1;

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .8;

% Model color
color = {'r','r'};
lw = .5;

% --------------------------------
% Load the ManualFit Parameters
load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);
PAR.image_filter = ManualFit.ImageData.image_filter;
PAR.imagepath = ManualFit.ImageData.imagepath;

%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
%cisco: analyze'em all
%frames = ManualFit.ImageData.frames4analysis; % number of frames in movie
%exp035:
%frames = [7 171];

%exp031:
%frames = [22 142];

%exp082:
%frames = [25 300];

%exp083:
%frames = [380 564]

%exp093:
%frames = [14 160];

%exp95: 
%frames  = [430 460];

%exp098:
%frames = [325 536];

%exp099:
%frames = [90 446];

%exp101:
frames = [90 445];

%exp104:
%frames = [915 1016];

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
ttime = (0:length(movidx)-1)*PAR.dt*1000;

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;


soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);

if interptag == 1
    for i=1:length(movidx)
        load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(movidx(i)) ...
            '.mat']);
        soln1(i,:) = xh(1:PAR.numfly*PAR.statedim)';
        clear xh InternalVariablesDS
    end
end



% %Now, perform interpolation to calculate the state at in between frames
% for k = 1:size(soln1,2)
%     soln1(:,k) = interp1(movidx,soln1(:,k),frames(1):frames(2),'spline')';
% end

%Now, perform interpolation to calculate the state at the frames where the
%estimated solution is bad
% fixframes = 78:81; %(exp100)
% fixframes = 120:133; %(exp082)
% fixidx = find(movidx == fixframes(1)) : find(movidx == fixframes(end));
% %Interpolate the displacements
% for k = 1:3
%     soln1(fixidx,k) = interp1(movidx,soln1(:,k),fixframes,'spline')';
% end
% 
% %Now interpolate the quaternion rotations
% bidx = 4:7;
% lwidx = 8:11;
% rwidx = 12:15;
% 
% %Body
% Qtmp = [soln1(fixidx(1)-1,bidx) ; soln1(fixidx(end)+1,bidx)]';
% QtmpNew = slerpMAT(Qtmp,length(fixframes)+1);
% soln1(fixidx,bidx) = QtmpNew(:,2:end-1)';
% 
% %Wings
% Qtmp = [soln1(fixidx(1)-1,lwidx) ; soln1(fixidx(end)+1,lwidx)]';
% QtmpNew = slerpMAT(Qtmp,length(fixframes)+1);
% soln1(fixidx,lwidx) = QtmpNew(:,2:end-1)';
% 
% Qtmp = [soln1(fixidx(1)-1,rwidx) ; soln1(fixidx(end)+1,rwidx)]';
% QtmpNew = slerpMAT(Qtmp,length(fixframes)+1);
% soln1(fixidx,rwidx) = QtmpNew(:,2:end-1)';



%%

%--------------------------------------------
%calculate theta-dot
theta_vec = sqrt(sum(soln1(:,4:6).^2,2));
theta_dot = gradient(theta_vec,PAR.dt);
dt = PAR.dt;

M2 = moviein(size(soln1));
for i= 1:size(soln1,1)
    kk = frames(1)+(i-1);

    % ------------------------------------
    % Load images
    %Number of digit places for the total number of frames
    %
    %cisco: need to change the definition of digits
    digits = length(num2str(numframes));
    for cam=1:PAR.numcam
        %load the raw grayscale image
        input_filename = sprintf(['%scam%03d/%s%0' num2str(digits) 'd%s'], ...
            PAR.imagepath, cam, PAR.stub, kk,PAR.image_filter(2:end));

        Input(:,:,cam) = imread(input_filename);
    end
    
    PAR.imgres = size(Input(:,:,1));
    PAR.imgres_crop = PAR.imgres;

%     load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(kk) ...
%         '.mat']);
%     load([PAR.solutionpath 'Features_' PAR.stub '/Features' num2str(kk) ...
%         '.mat']);


    clear flymod
    [x,y,z] = flymodQ(soln1(i,:),PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    %----------------------
    % Plot Movie Frame or Rendered model image
    %----------------------
    pts = [];
    for j = 1:length(x)
        pts{j} = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
        
        if j == 1
            % get dorsal edge of the body
            dorsalpts = [x{j}(:,10),y{j}(:,10),z{j}(:,10)];

        end
    
    end

    % Initialize projected planar points
    u = cell(PAR.numcam,length(x));
    v = cell(PAR.numcam,length(x));
    uT = cell(PAR.numcam,length(x));
    vT = cell(PAR.numcam,length(x));
    
    %---------------------------------------------------------
    % Calculate the spatial velocity
    skew = @(p) [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];
    scr = soln1(i,1:6)';
    G = screw2homo(scr);
    TT = G(:,4);
    theta = norm(scr(4:6));
    om = scr(4:6)./theta;
    %This is spatial velocity of body frame
    V = theta_dot(i).*[skew(om) scr(1:3);zeros(1,4)] * TT;
    V = V(1:3);
    scal = norm(V);
    Vel_pts = [TT(1:3)' ; TT(1:3)'+ 0.01.*V'];
    
    
    
    
    figure(fignum);
    clf;

    im = reshape(Input,PAR.imgres(1),[]);
    imagesc(im);
    colormap gray
    ax1 = gca;
    set(ax1,'units','pixels','visible','off','position',[0 0 PAR.imgres(2)*PAR.numcam PAR.imgres(1)].*plotscale);
    axis image

    %aspectratio = get(ax1,'PlotBoxAspectRatio');
    impos = get(ax1,'position');
    %we want position to be for single image, not composite image
    impos(3) = impos(3)/PAR.numcam;
    for k = 1:PAR.numcam
        ax2 = axes('position',[0 0 .5 .5]);
        hold on

        DLTparams = PAR.DLT(:,k);
        
        Vel_2Dpts = dlt_3D_to_2D(DLTparams,Vel_pts);
        
        for j = 1:length(pts)
            dorsal2Dpts = dlt_3D_to_2D(DLTparams,dorsalpts);
            uv = dlt_3D_to_2D(DLTparams,pts{j});

            u{k,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2));
            v{k,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2));

            uT{k,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2))';
            vT{k,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2))';
            
            if j == 1
                %Plot dorsal edge
                plot(ax2,dorsal2Dpts(:,1),dorsal2Dpts(:,2),'g-','linewidth',5);

                hold on;

                %Plot Head and tail
                %plot(ax2,uv([1 end],1),uv([1 end],2),'rs');
            end

            hold on;
            patch(u{k,j},v{k,j},[1 1 1],'facecolor','none','edgecolor',color{1},'linewidth',lw);
            patch(uT{k,j},vT{k,j},[1 1 1],'facecolor','none','edgecolor',color{1},'linewidth',lw);
        end
        
        %plot(Vel_2Dpts(:,1),Vel_2Dpts(:,2),'r.-','markersize',10,'linewidth',2);
        %keyboard
        
        set(ax2,'units','pixels','fontsize',12,'position',impos+[PAR.imgres(2)*(k-1) 0 0 0]*plotscale,'color','none','xlim',...
            [0.5 PAR.imgres(2)+0.5],'ylim',[-0.5 PAR.imgres(1)-0.5],'visible','off',...
            'xdir','normal','ydir','reverse');
        
        %Plot frame # on the 2nd camera image (in the middle)
        if k == 2
%             text(PAR.imgres(2)/2 - 50,20,[ '\color{white}Frame ' num2str(kk)], ...
%                 'fontsize',20);
             text(PAR.imgres(2)/2 - 50,20,['\color{black}' num2str(ttime(i),'%-3.2f') ' \color{black} ms'], ...
                'fontsize',20);
        end
    end

    figure(fignum);
    set(fignum,'units','pixels','position',[10 150 PAR.imgres(2)*PAR.numcam PAR.imgres(1)].*plotscale);
    %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])
    

    M2(i) = getframe(fignum);
%     M2(kk) = getframe(fignum);
end

