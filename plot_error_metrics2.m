%plot_error_metrics.m
%
%plot_error_metrics plots the state variables of a specified sequence over a
%specified number of frames.
%It will plot the estimated values and/or values from a prior database.
%
clear all;
%close all;
clc;

PAR.videopath = '../FlyVideo/';

% cutoff frequency for butterworth filter used to smooth time traces
Fcutoff_wing = 1000;
Fcutoff_body = 250;
PAR.dt = 1/6000;
Fs = 1/PAR.dt;

PAR.filetag = {'exp098000','exp099000','exp100000','exp101000','exp102000','exp104000'};

PAR.solutionpath = [PAR.videopath '/solutions/'];
frgroup = {[325 519],[90 425],[40 170],[90 402],[61 390],[800 1016]};

%These are the frames with regular wing beats used to fit a line through 
strokeReg = {[325 519],[270 500],[40 170],[196 396],...
             [100 390],[800 1016],[380 545],[25 135],...
             [440 640],[60 100],[1,100]};
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
splineMe = 0;

% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;

for m = 1%:length(PAR.filetag)
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
    %PAR.stub = 'exp098';
    %load(['kine/SavedKinematics/from_gwyneth/' PAR.stub '_fixed_eb.mat']);
    load(['kine/SavedKinematics/from_gwyneth/exp098_fixed_eb.mat']);

    %State = zeros(PAR.statedim,length(frames1));

    BodyTrans = zeros(3,frames(2));
    BodyAng = zeros(3,frames(2));

    for k = 1:numframes

        %% Get the body orientation
        q_body = data.kine.body.data.quat(:,movidx(k));
        % Get Translation
        T_body = data.kine.body.data.v_trans(:,movidx(k));
if 0       
        if (all(q_body == zeros(4,1)) || all(T_body == zeros(3,1)))
            addpath 'kine/kine_math';
            %Kine didn't calculate the quaternion at this frame
            %fill it in
            alpha = (pi/180) * data.kine.body.data.params(1,movidx(k));
            alpha = 2*pi - alpha;
            
            % Get length, azimuth, elevation for digitized anchor points
            % Correct direction for vector is HEAD - TAIL
            c1 = data.kine.body.data.coords(1,:,movidx(k));
            c2 = data.kine.body.data.coords(2,:,movidx(k));
            dc = c1 - c2;
            [psi,theta,c_length] = cart2sph(dc(1),dc(2),dc(3));
            
            % Get correct angles to define LAB to BODY rotation
            alpha = alpha;%pi-alpha;%(2*pi)-alpha; %DEBUG NEW SYSTEM: for new will just be alpha
            theta = -theta; % reverse theta because cart2sph defines elevation above the xy-plane as positive, but this is a negative rotation around the y-axis
            psi = psi; % don't reverse because rotate from body to lab frame

            % Calculate quaternion, rotate, scale, translate
            q_body = eulzyx2quat(alpha,theta,psi);
            
            data.kine.body.data.quat(:,movidx(k)) = q_body;
            
            %Now fix Translation
            model_3d = data.kine.body.config.model_coords;
            
            % Get Model length (really this should always be 1)
            if isfield(data.kine.body.config,'model_length')
                m_length = data.kine.body.config.model_length;
            else
                anch_1 = data.kine.body.config.anchor_array{1,2};
                anch_2 = data.kine.body.config.anchor_array{2,2};

                m1 = data.kine.body.config.model_coords(:,anch_1);
                m2 = data.kine.body.config.model_coords(:,anch_2);
                dm = m2 - m1;

                [m_theta,m_phi,m_length] = cart2sph(dm(1),dm(2),dm(3));
            end
            % Scale factor is digitized units/model units
            s = c_length / m_length;
            m_fit = quatRot_bod2lab(q_body, model_3d);
            m_fit = scale3D(m_fit,s);

            pt_name = data.kine.body.config.anchor_array{1,1}; % use first anchor point
            pt_col = data.kine.body.config.anchor_array{1,2};
            kine_row = strmatch(pt_name, data.kine.body.config.points);
            vb_anch = m_fit(:,pt_col);
            vf_anch = data.kine.body.data.coords(kine_row,:,movidx(movidx(k)));

            [m_fit, v_trans] = trans3D(m_fit, vb_anch, vf_anch');
            
            data.kine.body.data.v_trans(:,movidx(k)) = v_trans;
            
            save(['kine/SavedKinematics/from_gwyneth/' PAR.stub '_fixed_eb.mat'],'data');
        end
end      
        %I store the scalar part of quaternion at the end
        q_body = [q_body(2:4); q_body(1)];

        % add a rotation by pi along the roll axis since our model has the body
        % fixed frame with z-axis pointing ventral.
        q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);

        % Get Translation
        T_body = data.kine.body.data.v_trans(:,movidx(k));
        %Calculation translation vector -q_body*T10*q_body` + T_body
        %T10 is the approximate location of tailpt in our generative model
        %It is the next to last cross-section of model at the dorsal edge
%%%%%%%%%
        %cisco: unnecessary here(?)
%%%%%%%%%        
        %[xbody,ybody,zbody,s,th,X,Frenet] = flybodymod(ManualFit.params.bodyctr,...
        %    ManualFit.params.bodyrad,ManualFit.params.bodylen,PAR);

        %Say T10 is right on the tip of the thorax
        %T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
        %T10 = -T10;
        %[T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),T_body,q_body);

        %Calculate the equivalent location using manually fit body length as the
        %T_body used in our model for proper error comparison.
        bodyvec = diff(data.kine.body.data.coords(:,:,movidx(k)),1,1);
        bodyvec = -bodyvec./norm(bodyvec);
        T_body = T_body + .4226*data.kine.body.data.length(movidx(k)).*bodyvec';


        BodyTrans(:,movidx(k)) = T_body;

        
        %% Calculate the twist representation for these quaternions
        R_body = quat2matNEW(q_body);
        BodyAng(:,movidx(k)) = Rot2Euler(R_body');     
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % auto_init.m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %===============================================================
        % Get the wing orientations
        %% Left Wing
        q_Lwing(:,k) = data.kine.left_wing.data.quat(:,movidx(k));

        %I store the scalar part of quaternion at the end
        q_Lwing(:,k) = [q_Lwing(2:4,k); q_Lwing(1,k)];

        % Premultiply rotation by this alignment quaternion that takes into account
        % the orientation of coordinate axis fixed to the left wing.
        q_Lwingaxisalign = quat2matNEW([0 -1 0;1 0 0;0 0 1]);
        q_Lwing(:,k) = quatprod(q_Lwing(:,k),q_Lwingaxisalign);


        T_L(:,k) = data.kine.left_wing.data.v_trans(:,movidx(k));

        % I will calculate the relative rotations
        % for the wings by multiplying the quaternions directly.
        % Just take the orientation because I assume that the wing is fixed at the
        % joint.
        q_Lwing_rel(:,k) = quatprod([-q_body(1:3) ; q_body(4)],q_Lwing(:,k));

        % Relative translation from body fixed axis to wing joint
        params.T_Lwing_rel(:,k) = qxform([-q_body(1:3) ; q_body(4)],T_L(:,k) - T_body);      
        
        %--------------------------------------------
        %% Right Wing
        q_Rwing(:,k) = data.kine.right_wing.data.quat(:,movidx(k));

        %I store the scalar part of quaternion at the end
        q_Rwing(:,k) = [q_Rwing(2:4,k); q_Rwing(1,k)];

        % Premultiply rotation by this alignment quaternion that takes into account
        % the orientation of coordinate axis fixed to the right wing.
        q_Rwingaxisalign = quat2matNEW([0 1 0;1 0 0;0 0 -1]);
        q_Rwing(:,k) = quatprod(q_Rwing(:,k),q_Rwingaxisalign);


        T_R(:,k) = data.kine.right_wing.data.v_trans(:,movidx(k));

        %Just take the orientation because I assume that the wing is fixed at the
        %joint.
        q_Rwing_rel(:,k) = quatprod([-q_body(1:3); q_body(4)],q_Rwing(:,k));

        % Relative translation from body fixed axis to wing joint
        params.T_Rwing_rel(:,k) = qxform([-q_body(1:3) ; q_body(4)],T_R(:,k) - T_body);       
        
        % If the fourth element in the wing quaternions is not
        % positive, then change the sign of the the wing quaternions 
        % to make it match with the database prior.  Remember, q and -q
        % map to the same point in SO(3).
        if q_Lwing_rel(4,k) < 0
            q_Lwing_rel(:,k) = -q_Lwing_rel(:,k);
        end

        if q_Rwing_rel(4,k) < 0
           q_Rwing_rel(:,k) = -q_Rwing_rel(:,k);
        end

        % My ordering for the twist axes corresponds to the Rxyz Euler angles.
        %
        Theta_Lwing(:,k) = Rot2Joint(quat2matNEW(q_Lwing_rel(:,k)));
        Theta_Rwing(:,k) = Rot2Joint(quat2matNEW(q_Rwing_rel(:,k)));
        
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   WingAnglesCalc.m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate the wing tip trajectory & leading edge orientation
        
        %==========================
        % Left Wing
        %quat = SOLN(12:15,:);
        quat = q_Lwing;
        pts = zeros(2,size(quat,2));
        LEpts = zeros(2,size(quat,2));

        for k = 1:size(pts,2)
            G = quat2matNEW(quat(:,k));
            % Scale the wing tip point
            scale = PAR.params.wingscale*PAR.params.wing_tip2tip;
            pts(:,k) = scale.*[G(1,2) ; G(3,2)];
        end

        % indices of the wing trajectory that are used to calculate the
        % regression
        regidx = strokeReg{m} - frgroup{m}(1) + 1;
        %regidx = strokeReg{m};


        Rpts = pts;
    
        %==========================
        % Right Wing
        %quat = SOLN(8:11,:);
        quat = q_Rwing;
        
        for k = 1:size(pts,2)
            G = quat2matNEW(quat(:,k));
            pts(:,k) = scale.*[-G(1,2) ; -G(3,2)];
        end

        Lpts = pts;     
        
        %%=====================================
        %% Mean Stroke Plane
        %
        %Only take the points from the wing tip trajectory where the wingbeats
        %are regular to fit the line to.
        xpts = [Rpts(1,regidx(1):regidx(2)) Lpts(1,regidx(1):regidx(2))];
        ypts = [Rpts(2,regidx(1):regidx(2)) Lpts(2,regidx(1):regidx(2))];

        %P = polyfit(xpts,ypts,1);

        V = princomp([xpts' ypts']);
        tvec = V(:,1);
        nvec = [-tvec(2) ; tvec(1)];

        %M = mean([xpts' ypts'],1);
        M = [0 0]; %Set to wing joint location
        P(1) = -nvec(1)/nvec(2);
        P(2) = ( nvec'*M' ) / nvec(2);

        %Calculate model points to plot model shadow
        % Shift model to wing joint centered frame
        BL = PAR.params.bodyscale*(PAR.params.bodylen+PAR.params.headlen);
        RJTrans = BL.*([0.2021 0.1055 -0.1477]);
        [x,y,z] = flymodQ([-RJTrans 0 0 0 1],PAR.params,PAR);

        %Take the points in the X-Z plane that correspond to the dorsal and
        %ventral edge.
        Bodypts = [ [x{1}(:,4) ; x{1}(:,10)] [z{1}(:,4) ; z{1}(:,10)] ];

        %Figure out the angle that the stroke plane makes with the x-axis
        vv = [1 P(1)];
        vv = vv./norm(vv);

        %========================
        % Apply an extra rotation about the fixed body axis to calculate the
        % wing rotations in the appropriate stroke plane
        % NomPitch =  acos(vv(1))*(180/pi);

        % Set the pitch fixed to conform with that measured in Charlie David paper for
        % hovering.
        NomPitch = 62;

        q_mod = [0; sin(0.5.*NomPitch.*(pi/180)) ; 0; cos(0.5*NomPitch.*(pi/180))];

        %%================================
        % Add this rotation about the y axis to the wing quaternions so that 
        % the mean stroke plane aligns with the horizontal.
        %Rquat = zeros(size(SOLN(12:15,:)));
        %Lquat = zeros(size(SOLN(8:11,:)));

        for k = 1:size(q_Rwing_rel,2)
            Rquat(:,k) = quatprod(q_mod,q_Rwing_rel(:,k));
            Lquat(:,k) = quatprod(q_mod,q_Lwing_rel(:,k));
        end


        %%======================
        % Calculate the new wing tip trajectory, leading edge orientation
        % and wing angles
        LEpts_R = zeros(2,size(Rquat,2));
        LEpts_L = zeros(2,size(Lquat,2));

        theta_R = zeros(1,size(pts,2));
        phi_R = theta_R;
        alpha_R = theta_R;

        theta_L = zeros(1,size(pts,2));
        phi_L = theta_L;
        alpha_L = theta_L;
        
        for k = 1:size(pts,2)
            %%==============
            % LEFT WING
            G = quat2matNEW(Lquat(:,k));

            if k < size(pts,2)
                G1 = quat2matNEW(Lquat(:,k+1));
            else
                G1 = quat2matNEW(Lquat(:,1));
            end
            % Scale the wing tip point
            scale = PAR.params.wingscale*PAR.params.wing_tip2tip;
            Lpts(:,k) = scale.*[-G(1,2) ; -G(3,2)];

            LEpts_L(:,k) = [G(1,1) ; G(3,1)];

            %calculate wing angles
            [phitmp,thetatmp,dum] = cart2sph(-G(1,2),-G(2,2),-G(3,2));

            %make 90 degrees dorsal and -90 ventral for stroke amplitude
            %make deviation positive for dorsal side
            phi_L(k) = -phitmp - pi/2;
            theta_L(k) = -thetatmp;


            R2(:,k) = -G(1:3,2)*scale;
            R3(:,k) = -G1(1:3,2)*scale;
            %vector that points from trailing edge to leading edge
            LE1(:,k) = G(1:3,1);
            Vel1(:,k) = R3(:,k) - R2(:,k);
            Vel1(:,k) = Vel1(:,k)./norm(Vel1(:,k));

            %Switch signs depending on upstroke or downstroke (i.e. look at x
            %component of Leading edge vector
            if Vel1(1,k) > 0
                alpha_L(k) = acos(LE1(:,k)'*Vel1(:,k));
            else
                alpha_L(k) = -acos(LE1(:,k)'*Vel1(:,k));
            end



            %%==============
            % RIGHT WING

            G = quat2matNEW(Rquat(:,k));
            if k < size(pts,2)
                G1 = quat2matNEW(Rquat(:,k+1));
            else
                G1 = quat2matNEW(Rquat(:,1));
            end
            % Scale the wing tip point
            scale = PAR.params.wingscale*PAR.params.wing_tip2tip;
            Rpts(:,k) = scale.*[G(1,2) ; G(3,2)];

            % .96 is chord distance in mm of wing model
            chordscale = PAR.params.wingscale*.96;
            LEpts_R(:,k) = [G(1,1) ; G(3,1)];

            %calculate wing angles?
            [phitmp,thetatmp,dum] = cart2sph(G(1,2),G(2,2),G(3,2));

            %make 90 degrees dorsal and -90 ventral for stroke amplitude
            %make deviation positive for dorsal side
            phi_R(k) = phitmp - pi/2;
            theta_R(k) = -thetatmp;


            R(:,k) = G(1:3,2)*scale;
            R1(:,k) = G1(1:3,2)*scale;
            %vector that points from trailing edge to leading edge
            LE(:,k) = G(1:3,1);
            Vel(:,k) = R1(:,k) - R(:,k);
            Vel(:,k) = Vel(:,k)./norm(Vel(:,k));

            %Switch signs depending on upstroke or downstroke (i.e. look at x
            %component of Leading edge vector
            if Vel(1,k) > 0
                alpha_R(k) = acos(LE(:,k)'*Vel(:,k));
            else
                alpha_R(k) = -acos(LE(:,k)'*Vel(:,k));
            end


        end

        %Unwrap amplitude discontinuities
        ii = find(phi_L < -pi);
        phi_L(ii) = phi_L(ii) + 2*pi;

        ii = find(phi_R < -pi);
        phi_R(ii) = phi_R(ii) + 2*pi;


if ~splineMe
        %================================
        % Now, smooth the wing and body angles
        % Cutoff frequency from fft on a sample signal
        filter_order = 4;
        [b a] = butter( filter_order,Fcutoff_body*(2/Fs));

        for k = 1:size(BodyAng_auto,1)
            BodyAng_auto(k,:) = filtfilt(b,a,BodyAng_auto(k,:));
            BodyAng(k,:) = filtfilt(b,a,BodyAng(k,:));

        end

        [b a] = butter( filter_order,Fcutoff_wing*(2/Fs));
        
        phi_R = filtfilt(b,a,phi_R);
        phi_L = filtfilt(b,a,phi_L);

        theta_R = filtfilt(b,a,theta_R);
        theta_L = filtfilt(b,a,theta_L);

        alpha_R = filtfilt(b,a,alpha_R);
        alpha_L = filtfilt(b,a,alpha_L);
end

    
%% Plotting 

    
%%    %===================
    %Plot wing angles
    %t_end = movidx(end);

    %===================
    %Plot Stroke Amplitude   
    
    tt = (0:length(movidx)-1)*PAR.dt*1000;
    %tt = movidx;
    t_end = round(tt(end));
if 1    
    SAMPR = figure; 
    hold on
    plot(tt,phi_R*(180/pi),'r-','linewidth',2);
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-90 90],'ytick',-90:45:90,'yticklabel',{'-90','','0','','90'});
    box off

    title('S.Amp_R');    

    SAMPL = figure; 
    hold on
    plot(tt,phi_L*(180/pi),'b-','linewidth',2);
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-90 90],'ytick',-90:45:90,'yticklabel',{'-90','','0','','90'});
    box off

    title('S.Amp_L'); 
    
    %===================
    %Plot Stroke Plane Deviation (w.r.t. 62deg)    
    
    SPDR = figure; 
    plot(tt,theta_R*(180/pi),'r-','linewidth',2);
    hold on
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-40 40],'ytick',-40:20:40,'yticklabel',{'-40','','0','','40'});
    box off
    
    title('S.P.D_R');    
    
    SPDL = figure; 
    plot(tt,theta_L*(180/pi),'b-','linewidth',2);
    hold on    
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-40 40],'ytick',-40:20:40,'yticklabel',{'-40','','0','','40'});
    box off
    
    title('S.P.D_L'); 
end    
    %===================
    %Plot Angle of Attack (geometric)
    
    AoAR = figure; 
    plot(tt(1:end-1),alpha_R(1:end-1)*(180/pi),'r-','linewidth',2);
    hold on
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-180 180],'ytick',-180:90:180,'yticklabel',{'-180','','0','','180'});
    box off
    title('A-o-A_R');    

    AoAL = figure; 
    plot(tt(1:end-1),alpha_L(1:end-1)*(180/pi),'b-','linewidth',2);
    hold on
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-180 180],'ytick',-180:90:180,'yticklabel',{'-180','','0','','180'});
    box off
    title('A-o-A_L');   
    
    AvgSampR(1,:) = phi_R;
    AvgSampL(1,:) = phi_L;
    AvgSPDR(1,:) = theta_R;
    AvgSPDL(1,:) = theta_L;
    AvgAoAR(1,:) = alpha_R;
    AvgAoAL(1,:) = alpha_L;
    
 if 0   
    figure; 
    plot(BodyAng_auto(1,:).*(180/pi),'m');hold on;
    plot(BodyAng_auto(2,:).*(180/pi),'c');
    plot(BodyAng_auto(3,:).*(180/pi),'g');
    plot(BodyAng(1,:).*(180/pi),'r');
    plot(BodyAng(2,:).*(180/pi),'b');
    plot(BodyAng(3,:).*(180/pi),'k');
    clear phi_R phi_L;
 end    
    if 1
    [phi_Rn,phi_Ln,theta_Rn,theta_Ln,alpha_Rn,alpha_Ln] = WingAnglesCalc2;
    phi_R = phi_Rn ;
    phi_L = phi_Ln ;
    theta_R = (theta_Rn + theta_R)/2;
    theta_L = (theta_Ln + theta_L)/2;
    alpha_R = (alpha_Rn + alpha_R)/2;
    alpha_L = (alpha_Ln + alpha_L)/2;
    
    AvgSampR(2,:) = phi_R;
    AvgSampL(2,:) = phi_L;
    AvgSPDR(2,:) = theta_R;
    AvgSPDL(2,:) = theta_L;
    AvgAoAR(2,:) = alpha_R;
    AvgAoAL(2,:) = alpha_L;
    
    meanPhiR = mean(AvgSampR).*180/pi;
    meanPhiL = mean(AvgSampL).*180/pi;
    stdPhiR = std(AvgSampR).*180/pi;
    stdPhiL = std(AvgSampL).*180/pi;
    meanThetaR = mean(AvgSPDR).*180/pi;
    meanThetaL = mean(AvgSPDL).*180/pi;
    stdThetaR = std(AvgSPDR).*180/pi;
    stdThetaL = std(AvgSPDL).*180/pi;
    meanAlphaR = mean(AvgAoAR).*180/pi;
    meanAlphaL = mean(AvgAoAL).*180/pi;
    stdAlphaR = std(AvgAoAR).*180/pi;
    stdAlphaL = std(AvgAoAL).*180/pi;
    
    AvgSamp = [AvgSampR;AvgSampL];
    AvgSPD = [AvgSPDR;AvgSPDL];
    AvgAoA = [AvgAoAR;AvgAoAL]; 
    
    ErrSamp = [abs(AvgSamp(2,:)-AvgSamp(1,:));
               abs(AvgSamp(4,:)-AvgSamp(3,:))];
    ErrSPD  = [abs(AvgSPD(2,:)-AvgSPD(1,:));
               abs(AvgSPD(4,:)-AvgSPD(3,:))];
    ErrAoA  = [abs(AvgAoA(2,:)-AvgAoA(1,:));
               abs(AvgAoA(4,:)-AvgAoA(3,:))];           

    
    AvgErrSampR = sqrt( mean(ErrSamp(1,:).^2) ).*180/pi    ;
    AvgErrSPDR = sqrt( mean(ErrSPD(1,:).^2) ).*180/pi;
    AvgErrAoAR = sqrt( mean(ErrAoA(1,:).^2) ).*180/pi;
    AvgErrSampL = sqrt( mean(ErrSamp(2,:).^2) ).*180/pi;
    AvgErrSPDL = sqrt( mean(ErrSPD(2,:).^2) ).*180/pi;
    AvgErrAoAL = sqrt( mean(ErrAoA(2,:).^2) ).*180/pi;
    stdErrSampR = std(ErrSamp(1,:)).*180/pi;
    stdErrSPDR = std(ErrSPD(1,:)).*180/pi;
    stdErrAoAR = std(ErrAoA(1,:)).*180/pi;
    stdErrSampL = std(ErrSamp(2,:)).*180/pi;
    stdErrSPDL = std(ErrSPD(2,:)).*180/pi;
    stdErrAoAL = std(ErrAoA(2,:)).*180/pi;
    
    ErrSamp = mean(ErrSamp);
    ErrSPD = mean(ErrSPD);
    ErrAoA = mean(ErrAoA);  

    

    
    AvgErrSamp = sqrt( mean(ErrSamp).^2 ).*180/pi
    AvgErrSPD = sqrt( mean(ErrSPD).^2 ).*180/pi
    AvgErrAoA = sqrt( mean(ErrAoA).^2 ).*180/pi  
    stdErrSamp = std(ErrSamp).*180/pi;
    stdErrSPD = std(ErrSPD).*180/pi;
    stdErrAoA = std(ErrAoA).*180/pi ;
    
    meanPhi = mean(AvgSamp).*180/pi;
    meanTheta = mean(AvgSPD).*180/pi;    
    meanAlpha = mean(AvgAoA).*180/pi;    
    stdPhi = std(AvgSamp).*180/pi;
    stdTheta = std(AvgSPD).*180/pi;    
    stdAlpha = std(AvgAoA).*180/pi; 
    
    %keyboard
    %%    %===================
    %Plot wing angles
    %t_end = movidx(end);
if 1
    %===================
    %Plot Stroke Amplitude   
    figure(SAMPR);     hold on;
    plot(tt,phi_R*(180/pi),'m--','linewidth',2);
    
    figure(SAMPL);     hold on;
    plot(tt,phi_L*(180/pi),'c--','linewidth',2);
    
    figure(SPDR);     hold on;
    plot(tt,theta_R*(180/pi),'m--','linewidth',2);

    figure(SPDL);     hold on;
    plot(tt,theta_L*(180/pi),'c--','linewidth',2);    
end    
    figure(AoAR);     hold on;
    plot(tt,alpha_R*(180/pi),'m--','linewidth',2);
    
    figure(AoAL);     hold on;
    plot(tt,alpha_L*(180/pi),'c--','linewidth',2);
    end
end