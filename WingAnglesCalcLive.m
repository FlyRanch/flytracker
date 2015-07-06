close all; clear all; clc

ij = 9;
checkitout = 1; %1 AMP
                %2 SPD
                %3 AoA

expname = ['WingAngles_exp098c';
           'WingAngles_exp101c';
           'WingAngles_exp099c';
           'WingAngles_exp100c';
           'WingAngles_exp102c';
           'WingAngles_exp104c';
           'WingAngles_exp083c';
           'WingAngles_exp035c';
           'WingAngles_exp095c';
           ];
expnum  = [36;
           39;
           37;
           38;
           40;
           42;
           24;
            6;
           33;
          ];
          
color = ['r:';'b:';'g:';'m:';'c:';'r:';'b:';'g:'];
                   
load gmc/Combo2Results_n.mat
load(expname(ij,:))
      
scrsz = get(0,'ScreenSize');


PAR.videopath = '../FlyVideo/';
switch ij
    case 1
        %exp098: Tested
        PAR.filetag = 'exp098000';
        frames = [325 536];
    case 2
        %exp101: Tested
        PAR.filetag = 'exp101000';
        frames = [150 396];
    case 3
        %exp099: Tested
        PAR.filetag = 'exp099000';
        frames = [270 500];
    case 4
        %exp100: Tested
        PAR.filetag = 'exp100000';
        frames = [40 170];       
    case 5
        %exp102: Tested
        PAR.filetag = 'exp102000';
        frames = [61 390];
    case 6
        %exp104: NWY
        PAR.filetag = 'exp104000';
        frames = [800 1016];        
    case 7
        %exp083: Intriguing
        PAR.filetag = 'exp083000';
        frames = [380 545];
    case 8
        %exp035: Tested
        PAR.filetag = 'exp035000';
        frames = [25 135];
    case 9
        %exp095: Tested
        PAR.filetag = 'exp095000';
        frames = [440 639];        
end

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


%Initialize figures
%

%SPD  = figure('Position',[1 1 scrsz(3)/6 scrsz(4)*0.4]);
BO   = figure(12);
set(gcf,'Position',[scrsz(3)/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.45]);     
%AoA  = figure('Position',[scrsz(3)/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.4]); 
%WBF  = figure('Position',[scrsz(3)*2/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.4]); 
%Amp  = figure('Position',[scrsz(3)*2/6 1  scrsz(3)/6 scrsz(4)*0.4]);
%GMC  = figure('Position',[scrsz(3)*3/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.4]); 

%%  %===================
%Plot wing angles
%t_end = movidx(end);

    
    if checkitout == 1
        %===================
        %Plot Stroke Amplitude   
        %tt = movidx;
        SAmp = figure(11);
        set(gcf,'Position',[1 scrsz(4)*0.6 scrsz(3)/6 scrsz(4)*0.45]); 
        
        t_end = round(tt(end));
        figure(SAmp); 
        subplot(4,1,1);plot(tt,phi_R*(180/pi),'r-','linewidth',2);
        hold on
        subplot(4,1,1);plot(tt,phi_L*(180/pi),'b-','linewidth',2);
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-110 110],'ytick',-110:55:110,'yticklabel',{'-110','','0','','110'});
        hold off
        
        
    elseif checkitout == 2
        %===================
        %Plot Stroke Plane Deviation (w.r.t. 62deg)    
        SPD = figure(11);
        set(gcf,'Position',[1 scrsz(4)*0.6 scrsz(3)/6 scrsz(4)*0.45]); 
        
        figure(SPD); 
        title('SPD');box off
        subplot(4,1,1);plot(tt,theta_R*(180/pi),'r-','linewidth',2);
        hold on
        subplot(4,1,1);plot(tt,theta_L*(180/pi),'b-','linewidth',2);
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-70 70],'ytick',-70:35:70,'yticklabel',{'-70','','0','','70'});
        hold off

    elseif checkitout == 3
        %===================
        %Plot Angle of Attack
        AoA = figure(11);
        set(gcf,'Position',[1 scrsz(4)*0.6 scrsz(3)/6 scrsz(4)*0.45]); 
        
        figure(AoA)
        subplot(4,1,1);plot(tt(1:end-1),alpha_R(1:end-1)*(180/pi),'r-','linewidth',2);
        hold on
        subplot(4,1,1);plot(tt(1:end-1),alpha_L(1:end-1)*(180/pi),'b-','linewidth',2);
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-180 180],'ytick',-180:90:180,'yticklabel',{'-180','','0','','180'});
        hold off
    end
    
%===================
%Plot Body orientation
figure(BO);
plot(tt,BodyAng_auto(1,:)*(180/pi),'m-','linewidth',1);
hold on;
plot(tt,BodyAng_auto(2,:)*(180/pi),'y-','linewidth',1);
plot(tt,(BodyAng_auto(3,:) - BodyAng_auto(3,1))*(180/pi),'c-','linewidth',1);
set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
box off

figure(BO);     
title('B.O:E - roll(m), pitch(y), yaw(c)');box off

%===================
%Plot'em wingstrokes
[maxAmpR, minAmpR] = peakdet(phi_R, 0.1, tt);
[maxAmpL, minAmpL] = peakdet(phi_L, 0.1, tt);

figure(11); 
title('Stroke Amplitude');box off
subplot(4,1,4);hold on
    line([minAmpR(:,1)';minAmpR(:,1)'], [-180.*ones(1,length(minAmpR(:,1)));180.*ones(1,length(minAmpR(:,1)))],'Color','m','LineWidth',1);
subplot(4,1,4);hold on
    line([minAmpL(:,1)';minAmpL(:,1)'], [-180.*ones(1,length(minAmpL(:,1)));180.*ones(1,length(minAmpL(:,1)))],'Color','k','LineWidth',1);
%subplot(4,1,4);hold on
%    line([maxAmpR(:,1)';maxAmpR(:,1)'], [-180.*ones(1,length(maxAmpR(:,1)));180.*ones(1,length(maxAmpR(:,1)))],'Color','c','LineWidth',1);
%subplot(4,1,4);hold on
%    line([maxAmpL(:,1)';maxAmpL(:,1)'], [-180.*ones(1,length(maxAmpL(:,1)));180.*ones(1,length(maxAmpL(:,1)))],'Color','y','LineWidth',1);


%--------------------------------------------
%calculate theta-dot
theta_vec = sqrt(sum(soln1(:,4:6).^2,2));
theta_dot = gradient(theta_vec,PAR.dt);
dt = PAR.dt;

%M2 = moviein(size(soln1));

shazam = 0;
samplerate = 1;
movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);

for i= 1:size(soln1,1)
    kk = frames(1)+(i-1);

    % ------------------------------------
    % Load images
    %Number of digit places for the total number of frames
    %
    %cisco: need to change the definition of digits
    if (kk >999) && (shazam ~=1)
       PAR.stub = PAR.stub(1:length(PAR.stub)-1);
       shazam = 1;
    end
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
    set(fignum,'units','pixels','position',[1 1 PAR.imgres(2)*PAR.numcam PAR.imgres(1)].*plotscale);
    %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])
    
    
    %M2(i) = getframe(fignum);
    
    if checkitout == 1

        diffAngle = 3;

        ampDiff = (phi_R(i)-phi_L(i))*180/pi;
        if ampDiff >= diffAngle
            ampDiffColor = 'rx';
        elseif ampDiff <= -diffAngle
            ampDiffColor = 'bx';
        else
            ampDiffColor = 'kx';
        end

        figure(SAmp);
        subplot(4,1,2);plot(ttime(i),phi_R(i)*180/pi,'rx');hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-110 110],'ytick',-110:55:110,'yticklabel',{'-110','','0','','110'});    
        subplot(4,1,3);plot(ttime(i),phi_L(i)*180/pi,'bx');hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-110 110],'ytick',-110:55:110,'yticklabel',{'-110','','0','','110'});    
        subplot(4,1,4);plot(ttime(i),ampDiff,ampDiffColor);hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-20 20],'ytick',-20:10:20,'yticklabel',{'-20','','0','','20'});    

        figure(BO);
        plot(tt(i),BodyAng_auto(1,i)*(180/pi),ampDiffColor);

    %     M2(kk) = getframe(fignum);
    elseif checkitout == 2

        diffAngle = 0;

        devDiff = (phi_R(i)-phi_L(i))*180/pi;
        if devDiff >= diffAngle
            devDiffColor = 'rx';
        elseif devDiff <= -diffAngle
            devDiffColor = 'bx';
        else
            devDiffColor = 'kx';
        end

        figure(SPD);
        subplot(4,1,2);plot(ttime(i),theta_R(i)*180/pi,'rx');hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-70 70],'ytick',-70:35:70,'yticklabel',{'-70','','0','','70'}); 
        subplot(4,1,3);plot(ttime(i),theta_L(i)*180/pi,'bx');hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-70 70],'ytick',-70:35:70,'yticklabel',{'-70','','0','','70'});  
        subplot(4,1,4);plot(ttime(i),devDiff,devDiffColor);hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-50 50],'ytick',-50:25:50,'yticklabel',{'-50','','0','','50'});    

        figure(BO);
        %roll???
        %plot(tt(i),BodyAng_auto(1,i)*(180/pi),ampDiffColor);
        %yaw???
        %plot(tt(i),(BodyAng_auto(3,i) - BodyAng_auto(3,1))*(180/pi),devDiffColor);
        %pitch???
        plot(tt(i),BodyAng_auto(2,i)*(180/pi),devDiffColor);       
        
    elseif checkitout == 3

        diffAngle = 0;

        aoaDiff = (alpha_R(i)-alpha_L(i))*180/pi;
        if aoaDiff >= diffAngle
            aoaDiffColor = 'rx';
        elseif aoaDiff <= -diffAngle
            aoaDiffColor = 'bx';
        else
            aoaDiffColor = 'kx';
        end    
        
        figure(AoA)
        subplot(4,1,2);plot(ttime(i),alpha_R(i)*180/pi,'rx');hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-180 180],'ytick',-180:90:180,'yticklabel',{'-180','','0','','180'}); 
        subplot(4,1,3);plot(ttime(i),alpha_L(i)*180/pi,'bx');hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-180 180],'ytick',-180:90:180,'yticklabel',{'-180','','0','','180'});  
        subplot(4,1,4);plot(ttime(i),aoaDiff,aoaDiffColor);hold on
        set(gca,'plotboxaspectratio',[5 1 1],'xlim',[tt(1) t_end]);
        set(gca,'ylim',[-100 100],'ytick',-100:50:100,'yticklabel',{'-100','','0','','100'}); 
        
        figure(BO);
        %roll???
        %plot(tt(i),BodyAng_auto(1,i)*(180/pi),ampDiffColor);
        %yaw???
        %plot(tt(i),(BodyAng_auto(3,i) - BodyAng_auto(3,1))*(180/pi),devDiffColor);
        %pitch???
        plot(tt(i),BodyAng_auto(2,i)*(180/pi),aoaDiffColor);       
                
        
    end
end
