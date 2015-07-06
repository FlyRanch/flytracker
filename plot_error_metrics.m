%plot_error_metrics.m
%
%plot_error_metrics plots the state variables of a specified sequence over a
%specified number of frames.
%It will plot the estimated values and/or values from a prior database.
%
clear all
close all

PAR.videopath = 'video/';

PAR.filetag = {'exp098','exp099','exp100','exp101','exp102','exp104'};

PAR.solutionpath = [PAR.videopath '/solutions/'];
frgroup = {[325 521],[90 425],[40 170],[90 402],[61 390],[800 1016]};
%frgroup = {[325 521],[90 425],[40 170],[90 402],[800 1016]};

%This is the start frame for indexing every five points except exp102 (3
%pts)
FrameStart = [326 96 41 97 70 802];
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

for m = 4%1:length(PAR.filetag)
    PAR.stub = [PAR.filetag{m}];
    fignum = 1;

    % --------------------------------
    % Load the ManualFit Parameters
    load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);

    %----------------------------------------------------
    %load video data into buffer
    %----------------------------------------------------
    frames = frgroup{m}; % number of frames in movie

    samplerate = 1;
    movidx = frames(1):samplerate:frames(2);
    numframes = length(movidx);

    solidx = [1 frames(2)];
    
    sampidx = FrameStart(m):5:frames(2);
    if sampidx(1) ~= frames(1)
        sampidx = [frames(1) sampidx];
    end
    if sampidx(end) ~= frames(2)
        sampidx = [sampidx frames(2)];
    end
    
    % Assign model parameters
    PAR.params = ManualFit.params;
    PAR.DLT = ManualFit.DLT;
    PAR.cam = ManualFit.cam;


    soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
    SOLN = zeros(length(solidx(1):solidx(2)),PAR.numfly*PAR.statedim);
    
    if m <0%== 5
        soln1 = zeros(length(movidx),PAR.numfly*12);
        SOLN = zeros(length(solidx(1):solidx(2)),PAR.numfly*12);
    end
    
    if interptag == 1
        for i=1:length(movidx)
            load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(movidx(i)) ...
                '.mat']);
            if m <0%== 5
                soln1(i,:) = xh(1:PAR.numfly*12)';
            else
                soln1(i,:) = xh(1:PAR.numfly*PAR.statedim)';
            end
            clear xh InternalVariablesDS
        end
    end

    %SOLN = soln1;

    %Now, perform interpolation to calculate the state at in between frames
    for k = 1:size(soln1,2)
        SOLN(frames(1):frames(2),k) = interp1(movidx,soln1(:,k),frames(1):frames(2),'spline')';
    end
    
    SOLN = SOLN';
    
    
%%  ===============================================================
    % Okay, now I need to calculate the points associated with each state.
    % Take a subset of those points, and spline them, then, recalculate
    % the body quaternion 
    subSOLN = SOLN(1:7,sampidx);
    clear x y z xx yy zz T
    for k = 1:size(subSOLN,2)
        [X,Y,Z] = flymodQ(subSOLN(:,k),PAR.params,PAR);
        % Get pts
        x(:,k) = [X{1}(end,4) %head x-axis
            X{1}(8,7)%lateral y-axis
            X{1}(8,4)];%ventral z-axis
        y(:,k) = [Y{1}(end,4)
            Y{1}(8,7)
            Y{1}(8,4)];
        z(:,k) = [Z{1}(end,4)
            Z{1}(8,7)
            Z{1}(8,4)];
    end
    
    %spline the points
    for k = 1:size(x,1)
        xx(k,:) = spline(sampidx,x(k,:),movidx);
        yy(k,:) = spline(sampidx,y(k,:),movidx);
        zz(k,:) = spline(sampidx,z(k,:),movidx);
    end
    
    %spline the Translation
    T(1,:) = spline(sampidx,SOLN(1,sampidx),movidx);
    T(2,:) = spline(sampidx,SOLN(2,sampidx),movidx);
    T(3,:) = spline(sampidx,SOLN(3,sampidx),movidx);
    
    xv = [xx(1,:) ; yy(1,:) ; zz(1,:)] - T;
    mag = sqrt(sum(xv.^2,1));
    mag = repmat(mag,3,1);
    Rx = xv./mag;
    
    yv = [xx(2,:) ; yy(2,:) ; zz(2,:)] - T;
    mag = sqrt(sum(yv.^2,1));
    mag = repmat(mag,3,1);
    Ry = yv./mag;
    
    zv = [xx(3,:) ; yy(3,:) ; zz(3,:)] - T;
    mag = sqrt(sum(zv.^2,1));
    mag = repmat(mag,3,1);
    Rz = zv./mag;
    
    
%% ===============================================================

    % Convert the automatically captured body data into Euler components
    BodyTrans_auto = zeros(3,frames(2));
    BodyAng_auto = zeros(3,frames(2));
    
    if m <0%== 5
        for k = frames(1):frames(2)
            G = screw2homo(SOLN(1:6,k));
            BodyAng_auto(:,k) = Rot2Euler(G(1:3,1:3)');
            BodyTrans_auto(:,k) = G(1:3,4);
        end
    else
        for k = frames(1):frames(2)
            %G = screw2homo(SOLN(k,1:6));
            %Rbod = quat2matNEW(SOLN(4:7,k));
            Rbod = [Rx(:,k-frames(1)+1) Ry(:,k-frames(1)+1) Rz(:,k-frames(1)+1)];
            BodyAng_auto(:,k) = Rot2Euler(Rbod');%Rot2Euler(G(1:3,1:3)');
            BodyTrans_auto(:,k) = T(:,k-frames(1)+1);%SOLN(1:3,k);%G(1:3,4);
        end
    end

    

    %% Load the Manually Tracked data
    load(['kine/SavedKinematics/from_gwyneth/' PAR.stub '_fixed_eb.mat']);

    %State = zeros(PAR.statedim,length(frames1));

    BodyTrans = zeros(3,frames(2));
    BodyAng = zeros(3,frames(2));

    for k = movidx

        %% Get the body orientation
        q_body = data.kine.body.data.quat(:,k);
        % Get Translation
        T_body = data.kine.body.data.v_trans(:,k);
        
%         if (all(q_body == zeros(4,1)) || all(T_body == zeros(3,1)))
%             addpath 'kine/kine_math';
%             %Kine didn't calculate the quaternion at this frame
%             %fill it in
%             alpha = (pi/180) * data.kine.body.data.params(1,k);
%             alpha = 2*pi - alpha;
%             
%             % Get length, azimuth, elevation for digitized anchor points
%             % Correct direction for vector is HEAD - TAIL
%             c1 = data.kine.body.data.coords(1,:,k);
%             c2 = data.kine.body.data.coords(2,:,k);
%             dc = c1 - c2;
%             [psi,theta,c_length] = cart2sph(dc(1),dc(2),dc(3));
%             
%             % Get correct angles to define LAB to BODY rotation
%             alpha = alpha;%pi-alpha;%(2*pi)-alpha; %DEBUG NEW SYSTEM: for new will just be alpha
%             theta = -theta; % reverse theta because cart2sph defines elevation above the xy-plane as positive, but this is a negative rotation around the y-axis
%             psi = psi; % don't reverse because rotate from body to lab frame
% 
%             % Calculate quaternion, rotate, scale, translate
%             q_body = eulzyx2quat(alpha,theta,psi);
%             
%             data.kine.body.data.quat(:,k) = q_body;
%             
%             %Now fix Translation
%             model_3d = data.kine.body.config.model_coords;
%             
%             % Get Model length (really this should always be 1)
%             if isfield(data.kine.body.config,'model_length')
%                 m_length = data.kine.body.config.model_length;
%             else
%                 anch_1 = data.kine.body.config.anchor_array{1,2};
%                 anch_2 = data.kine.body.config.anchor_array{2,2};
% 
%                 m1 = data.kine.body.config.model_coords(:,anch_1);
%                 m2 = data.kine.body.config.model_coords(:,anch_2);
%                 dm = m2 - m1;
% 
%                 [m_theta,m_phi,m_length] = cart2sph(dm(1),dm(2),dm(3));
%             end
%             % Scale factor is digitized units/model units
%             s = c_length / m_length;
%             m_fit = quatRot_bod2lab(q_body, model_3d);
%             m_fit = scale3D(m_fit,s);
% 
%             pt_name = data.kine.body.config.anchor_array{1,1}; % use first anchor point
%             pt_col = data.kine.body.config.anchor_array{1,2};
%             kine_row = strmatch(pt_name, data.kine.body.config.points);
%             vb_anch = m_fit(:,pt_col);
%             vf_anch = data.kine.body.data.coords(kine_row,:,k);
% 
%             [m_fit, v_trans] = trans3D(m_fit, vb_anch, vf_anch');
%             
%             data.kine.body.data.v_trans(:,k) = v_trans;
%             
%             save(['kine/SavedKinematics/from_gwyneth/' PAR.stub '_fixed_eb.mat'],'data');
%         end
        
        %I store the scalar part of quaternion at the end
        q_body = [q_body(2:4); q_body(1)];

        % add a rotation by pi along the roll axis since our model has the body
        % fixed frame with z-axis pointing ventral.
        q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);

        % Get Translation
        T_body = data.kine.body.data.v_trans(:,k);
        %Calculation translation vector -q_body*T10*q_body` + T_body
        %T10 is the approximate location of tailpt in our generative model
        %It is the next to last cross-section of model at the dorsal edge
        [xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(ManualFit.params.bodyctr,...
            ManualFit.params.bodyrad,ManualFit.params.bodylen,PAR);

        %Say T10 is right on the tip of the thorax
        %T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
        %T10 = -T10;
        %[T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),T_body,q_body);

        %Calculate the equivalent location using manually fit body length as the
        %T_body used in our model for proper error comparison.
        bodyvec = diff(data.kine.body.data.coords(:,:,k),1,1);
        bodyvec = -bodyvec./norm(bodyvec);
        T_body = T_body + .4226*data.kine.body.data.length(k).*bodyvec';


        BodyTrans(:,k) = T_body;

        
        %% Calculate the twist representation for these quaternions
        R_body = quat2matNEW(q_body);
        BodyAng(:,k) = Rot2Euler(R_body');     
    end

    %Common frames between manual and automatically tracked data
    ff = frames(1):frames(2);
    BodyAng = BodyAng.*(180/pi);
    BodyAng_auto = BodyAng_auto.*(180/pi);
    tt = (0:length(ff)-1)*PAR.dt*1000;
    tsampidx = sampidx - frames(1) + 1;
    
    if m > 0
        figure;
        subplot(1,2,1);
        mkr = {'--','-','o','*'};
        hold on;
        
        % This is a hack to move the yaw direction so it fits on the graph better 
        BodyAng(3,:) = BodyAng(3,:) - 100;
        BodyAng_auto(3,:) = BodyAng_auto(3,:) - 100;
        
        %plot prior database
        hh = plot(tt,BodyAng(1,ff),['r' mkr{1}],tt,BodyAng(2,ff),['g' mkr{1}],tt,BodyAng(3,ff),['b' mkr{1}]);
        plot(tt(tsampidx),BodyAng(1,sampidx),['r' mkr{3}],tt(tsampidx),BodyAng(2,sampidx),['g' mkr{3}],tt(tsampidx),BodyAng(3,sampidx),['b' mkr{3}]);

        for i = 1:length(hh)
            set(hh(i),'linewidth',1)
        end

        
        %plot tracked data
        hh = plot(tt,BodyAng_auto(1,ff),['r' mkr{2}],tt,BodyAng_auto(2,ff),['g' mkr{2}],tt,BodyAng_auto(3,ff),['b' mkr{2}]);
        plot(tt(tsampidx),BodyAng_auto(1,sampidx),['r' mkr{4}],tt(tsampidx),BodyAng_auto(2,sampidx),['g' mkr{4}],tt(tsampidx),BodyAng_auto(3,sampidx),['b' mkr{4}]);

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
        plot(tt(tsampidx),BodyTrans(1,sampidx),['r' mkr{3}],tt(tsampidx),BodyTrans(2,sampidx),['g' mkr{3}],tt(tsampidx),BodyTrans(3,sampidx),['b' mkr{3}]);

        for i = 1:length(hh)
            set(hh(i),'linewidth',1)
        end
        
        
        %plot tracked data
        hh = plot(tt,BodyTrans_auto(1,ff),['r' mkr{2}],tt,BodyTrans_auto(2,ff),['g' mkr{2}],tt,BodyTrans_auto(3,ff),['b' mkr{2}]);
        plot(tt(tsampidx),BodyTrans_auto(1,sampidx),['r' mkr{4}],tt(tsampidx),BodyTrans_auto(2,sampidx),['g' mkr{4}],tt(tsampidx),BodyTrans_auto(3,sampidx),['b' mkr{4}]);

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
    end
    
    %Normalize by body length
    BodyTrans = BodyTrans ./ (ManualFit.params.bodyscale*(...
        ManualFit.params.bodylen + ManualFit.params.headlen));
    
    BodyTrans_auto = BodyTrans_auto ./ (ManualFit.params.bodyscale*(...
        ManualFit.params.bodylen + ManualFit.params.headlen));
    
    %% Error in body transformation
    AngErr{m} = (BodyAng_auto(:,ff) - BodyAng(:,ff));
    TransErr{m} = (BodyTrans_auto(:,ff) - BodyTrans(:,ff));
end

for k = 1:length(PAR.filetag)
    AngRMS(:,k) = sqrt(mean(AngErr{k}.^2,2));
    AngSTD(:,k) = std(AngErr{k},0,2);
    TransRMS(:,k) = sqrt(mean(TransErr{k}.^2,2));
    TransSTD(:,k) = std(TransErr{k},0,2);
    VecRMS(:,k) = mean( sqrt( sum(TransErr{k}.^2,1) ),2 );
    VecSTD(:,k) = std(sqrt( sum(TransErr{k}.^2,1) ),0,2);
end

figure; 
barweb(AngRMS,AngSTD,[],{'Roll','Pitch','Yaw'},'Body Angle Error',[],...
    'Degrees',[],'y',[]);

figure;
barweb([VecRMS],[VecSTD],[],[],'Body Translation Error',[],...
    'Body Length',[],'y',[]);