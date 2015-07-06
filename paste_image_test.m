% Paste_image_test.m
% 
% This file views the tracking solution at a particular frame
% This function tests the scaling of quaternion body motion

PAR.videopath = '../FlyVideo';
PAR.filetag = 'exp098';


PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.stub = [PAR.filetag];
fignum = 1;
frame = 100;
solidx = [frame:1:frame];
solidx(end) = frame;

for k = 1:3
%   close(10*fignum+k)
end

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .95;

plotpredict = 0;
plotcur = 1;
cptsTag = 0;   %Set to 1 to see the corresponding points
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
PAR.L3 = 30; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing
PAR.T3 = 3; %# of steps along thickness of wing

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
clear Input
for cam=1:PAR.numcam
    input_filename = sprintf('%s/%s/cam%03d/%s%06d.bmp', ...
        PAR.videopath,PAR.filetag,cam,PAR.filetag,frame);
    Input(:,:,cam) = imread(input_filename);
end

PAR.imgres = size(Input(:,:,1));
PAR.imgres_crop = PAR.imgres;

% ------------------------------------
% Plot images
for i = 1:PAR.numcam
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
end
%------------------------------------------
% Overlay the Model Mesh and other Features as described
%------------------------------------------

%Changing colors over time
%map = colormap('copper');
map = [0 0 0
    1 1 0];


load([PAR.solutionpath 'fly_' PAR.stub '_BodyPredictQ/fly' num2str(frame) '.mat']);
load([PAR.solutionpath 'Features_' PAR.stub '_BodyPredictQ/Features' num2str(frame) ...
    '.mat']);

scale = [.8:.05:1.1];
    
colors = [interp1(1:size(map,1),map(:,1),linspace(1,size(map,1),length(scale)))
    interp1(1:size(map,1),map(:,2),linspace(1,size(map,1),length(scale)))
    interp1(1:size(map,1),map(:,3),linspace(1,size(map,1),length(scale)))]';

for r = 1:length(scale)
    
    % ------------------------------------
    % Load the solution data
%     load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(frame) '.mat']);
%     load([PAR.solutionpath 'Features_' PAR.stub '/Features' num2str(frame) ...
%         '.mat']);

    Tb = xh(1:3);
    qb = xh(4:7);
    
    theta = 2*acos(qb(4));
    omega = qb(1:3)./sin(theta/2);
    
    beta = sin(theta);
    gamma = 1-cos(theta);
    omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        
    R = eye(3) + omegav*beta + omegav*omegav*gamma;
    A = (eye(3) - R)*omegav + omega*omega'.*theta;
    
    v = inv(A)*Tb;
    
    %%-------------------
    % scaled body motion
    ttheta = scale(r)*theta;
    
    qbnew = [omega.*sin(ttheta/2) ; cos(ttheta/2)];
    
    beta = sin(ttheta);
    gamma = 1-cos(ttheta);
        
    R = eye(3) + omegav*beta + omegav*omegav*gamma;
    
    %Now compute the Translation vector
    Tnew = ( (eye(3) - R)*omegav + omega*omega'.*ttheta )*v;
    
    sol = xh;
    sol(1:3) = Tnew;
    sol(4:7) = qbnew;


    clear flymodQ
    %clear flymod_BodyOnly
    [x,y,z] = flymodQ(sol,PAR.params,PAR);

    %----------------------
    % Plot Movie Frame or Rendered model image
    %----------------------
    pts = [];
    for i = 1:length(x)
        PAR.modsample(i) = size(x{i},1);
       

        pts{i} = [reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
        if i == 1
            % get dorsal edge of the body
            dorsalpts = [x{i}(:,10),y{i}(:,10),z{i}(:,10)];
        end

    end

    % Initialize projected planar points
    u = cell(PAR.numcam,length(x));
    v = cell(PAR.numcam,length(x));
    uT = cell(PAR.numcam,length(x));
    vT = cell(PAR.numcam,length(x));


    for i = 1:PAR.numcam

        DLTparams = PAR.DLT(:,i);

        shiftax = [-1 -1 -1];
        xax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 11 5 14]+ repmat(shiftax,2,1));
        yax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 6 14]+ repmat(shiftax,2,1));
        zax = dlt_3D_to_2D(DLTparams,[10 5 14 ; 10 5 15]+ repmat(shiftax,2,1));

        
        axes(subax(2*i));
        hold on;
        
        for j = 1:length(pts)

            dorsal2Dpts = dlt_3D_to_2D(DLTparams,dorsalpts);
            uv = dlt_3D_to_2D(DLTparams,pts{j});
            


            u{i,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2));
            v{i,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2));

            uT{i,j} = reshape(uv(:,1),size(x{j},1),size(x{j},2))';
            vT{i,j} = reshape(uv(:,2),size(x{j},1),size(x{j},2))';

            if j == 1
                %Plot dorsal edge
                %plot(gca,dorsal2Dpts(:,1),dorsal2Dpts(:,2),'r-','linewidth',5);

                hold on;

                %Plot Head and tail
                %plot(ax2,uv([1 end],1),uv([1 end],2),'rs');
            end


            if plotcur == 1
                %patch(u{i,j},v{i,j},[1 1 1],'facecolor','none','edgecolor','w');
                %patch(uT{i,j},vT{i,j},[1 1 1],'facecolor','none','edgecolor','w');
                patch(u{i,j},v{i,j},[1 1 1],'facecolor','none','edgecolor',colors(r,:));
                patch(uT{i,j},vT{i,j},[1 1 1],'facecolor','none','edgecolor',colors(r,:));
            end



        end

        %plot(xax(:,1),xax(:,2),'r.-',yax(:,1),yax(:,2),'g.-',zax(:,1),zax(:,2),'b.-')
        %label the axes
        %text(xax(2,1),xax(2,2),'X','fontsize',12)
        %text(yax(2,1),yax(2,2),'Y','fontsize',12)
        %text(zax(2,1),zax(2,2),'Z','fontsize',12)
        
        %--- Plot Correspondence between points
        if cptsTag == 1;
            %         datapts = Features(frame).DataptsFullIC{i};
            %         mdlpts = Features(frame).ModelptsIC{i};

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

        set(subax(2*i),'units','pixels','fontsize',12,'position',impos,'color','none','xlim',...
            [0.5 PAR.imgres(2)+0.5],'ylim',[-0.5 PAR.imgres(1)-0.5],'visible','off',...
            'xdir','normal','ydir','reverse');

        %set(fignum,'position',[0 500 PAR.imgres(2)*plotscale PAR.imgres(1)*plotscale])
        set(fignum*10+i,'units','pixels','position',[20+(i-1)*PAR.imgres(2) 50 PAR.imgres(2) PAR.imgres(1)].*plotscale);
        %title(['Frame ' num2str(frame) ', IC: p0 = ' num2str([sol(1:8)]) ])
        if r == length(solidx)
            text(PAR.imgres(2)/2 - 50,20,['\color{white}Frame ' num2str(frame)], ...
                'fontsize',20);
        end
    end
end

