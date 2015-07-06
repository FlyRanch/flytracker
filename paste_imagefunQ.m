function paste_imagefunQ(p,frame,ManualFit,PAR,varargin)
% Paste_image.m
% This file views the tracking solution at a particular frame


fignum = 1;


%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = 1;

plotpredict = 0;
plotcur = 1;
cptsTag = 0;   %Set to 1 to see the corresponding points
nrmlTag = 0;
renIMtag = 0; 



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

%----------------------------
% Evaluate Model at solution
%----------------------------
% Manual fit solution at first frame
% sol = reshape(ManualFit.soln,[],1);

% load testSigpts
% sol = testSigpts(:,10);

% Actual solution
sol = p;

sol1 = p;
    
clear flymodQ
%clear flymod_BodyOnly
[x,y,z] = flymodQ(sol,PAR.params,PAR);
[x1,y1,z1] = flymodQ(sol1,PAR.params,PAR);
for j = 1:length(x);
    PAR.modsample(j) = size(x{j},1);
end
%----------------------
% Plot Movie Frame or Rendered model image
%----------------------
pts = [];
for i = 1:length(x)
    pts{i} = [reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
    if i == 1
        % get dorsal edge of the body
        dorsalpts = [x{i}(:,10),y{i}(:,10),z{i}(:,10)];
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
            patch(u{i,j},v{i,j},[1 1 1],'facecolor','none','edgecolor','w');
            patch(uT{i,j},vT{i,j},[1 1 1],'facecolor','none','edgecolor','w');
        end
        if plotpredict == 1
            patch(u1{i,j},v1{i,j},[1 1 1],'facecolor','none','edgecolor','y');
            patch(uT1{i,j},vT1{i,j},[1 1 1],'facecolor','none','edgecolor','y');
        end
        
        
    end

    plot(xax(:,1),xax(:,2),'r.-',yax(:,1),yax(:,2),'g.-',zax(:,1),zax(:,2),'b.-')
    %label the axes
    text(xax(2,1),xax(2,2),'X','fontsize',12)
    text(yax(2,1),yax(2,2),'Y','fontsize',12)
    text(zax(2,1),zax(2,2),'Z','fontsize',12)

    %--- Plot Correspondence between points
    if cptsTag == 1;
        Features = varargin{1};
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
        keyboard
        [Y,idx,IMout,Nrml] = renderflyMOD2(pts,PAR,i);

        figure;
        imagesc(IMout); colormap gray;
        hold on;
        for j = 1:length(Y)
            plot(Y{j}(:,1),Y{j}(:,2),'b.');
            
            quiver(Y{j}(:,1),Y{j}(:,2),Nrml{j}(:,1),Nrml{j}(:,2),'r')
        end
    end
end