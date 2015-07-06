% Code to calculate the wing angles from the Fly Tracker results
%
% alpha - angle of attack
% phi - stroke deviation
% phi_dot - stroke deviation rate
%
% output -
% Wing Angles N x 3 x 2 matrix for N frames; Right wing is first element in 3rd dimension.
% Wing Angles(1,:,1) = [alpha phi phi_dot];
%
% Here I'll use the projection of the tip of the wing to calculate the
% stroke plane (From Altshuler, D. L.; et. al. "Short-amplitude high-frequency 
% wing strokes determine the aerodynamics of honeybee flight" PNAS, 2005.)
% By fitting an line through the wing trajectory

function [phi_R,phi_L,theta_R,theta_L,alpha_R,alpha_L] = WingAnglesCalc2

    PAR.videopath = '../FlyVideo/';
    PAR.filetag = {'exp098000','exp099000','exp100000','exp101000',...
                   'exp102000','exp104000','exp083000','exp035000',...
                   'exp095000','exp057000','exp031000'};

    PAR.solutionpath = [PAR.videopath '/solutions/'];
    %These are the frames where we calculate a tracker solution
    frgroup = {[325 519],[270 500],[40 170],[150 396],...
               [146 390],[800 1016],[380 545],[25 135],...
               [440 640],[1 120],[1,100]};

    %These are the frames with regular wing beats used to fit a line through 
    strokeReg = {[325 519],[270 500],[40 170],[196 396],...
                 [146 390],[800 1016],[380 545],[25 135],...
                 [440 640],[60 100],[1,100]};

    %exp101 -> [196 396]

    %exp102 07th wingbeat 300:324
    %       01st wingbeat 165:189
    %exp035 01st wingbeat 35:55
    %       02nd wingbeat 55:75
    %       05th wingbeat 115:134
    %exp083 01st wingbeat 432:450
    %       02nd wingbeat 450:473
    %       03rd wingbeat 473:496
    %       04th wingbeat 496:516
    %       05th wingbeat 516:537
    %exp098 01st wingbeat 340:368
    %       03rd wingbeat 394:421
    %       05th wingbeat 445:470 
    %       07th wingbeat 495:520
    %exp101 01st wingbeat 170:192
    %       05th wingbeat 263:287
    %       06th wingbeat 287:309
    %       10th wingbeat 375:395
    %       11th wingbeat 397:419
    %exp099 01st wingbeat 270:298
    %exp057 01st wingbeat 92:114
    %       02nd wingbeat 114:135
    %       03rd wingbeat 135:157

    %cisco: explanation?
    %exp100:  63:169;
    %exp102: 146:389;
    %exp104: 915:1015; %do-over?
    %exp083: 432:517;
    %exp035:  35:134;
    %exp095: 440:639;  %do-over?
    %exp098: 325:519;
    wingbeatframes = 325:519;

    % cutoff frequency for butterworth filter used to smooth time traces
    Fcutoff_wing = 1000;
    Fcutoff_body = 250;

    interptag = 1;

    PAR.pixpermm = 1;
    PAR.numfly = 1;

    %Number of parameters of the model (i.e. 8 control points)
    PAR.mdlpar = 15*ones(1,PAR.numfly);
    PAR.statedim = PAR.mdlpar;
    PAR.etamax = 0;

    %spline order
    PAR.c = 4;

    PAR.L1 = 40; %# of steps for body along length
    PAR.L2 = 15; %# of steps for head along length
    PAR.L3 = 20;%30; %# of steps for wing around the boundary
    PAR.T1 = 13; %# of theta steps for head and body
    PAR.T2 = 2; %# of steps towards center of wing
    PAR.T3 = 3; %# of steps along thickness of wing

    % - Camera Info
    PAR.dt = 1/6000;  %Framerate of the camera
    Fs = 1/PAR.dt;

    PAR.numcam = 3;

    % Set "m" equal to the index of PAR.filetag that corresponds to the video you want to analyze.
    for m = 1
        PAR.stub = [PAR.filetag{m}];

        % --------------------------------
        % Load the ManualFit Parameters
        load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub ]);

        %----------------------------------------------------
        %load solution data into buffer
        %----------------------------------------------------
        frames = frgroup{m}; % number of frames in movie
        samplerate = 1;
        movidx = frames(1):samplerate:frames(2);
        numframes = length(movidx);

        solidx = [1 frames(2)];


        % Assign model parameters
        PAR.params = ManualFit.params;
        PAR.DLT = ManualFit.DLT;
        PAR.cam = ManualFit.cam;

        soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
        SOLN = zeros(length(solidx(1):solidx(2)),PAR.numfly*PAR.statedim);

        if interptag == 1
            for i=1:length(movidx)
                load([PAR.solutionpath 'fly_' PAR.stub  '/fly' num2str(movidx(i)) ...
                    '.mat']);
                if m <0%== 5
                    soln1(i,:) = xh(1:PAR.numfly*12)';
                else
                    soln1(i,:) = xh(1:PAR.numfly*PAR.statedim)';
                end
                clear xh InternalVariablesDS
            end
        end

        SOLN = soln1;

        SOLN = SOLN';


        %%  ===============================================================
        % Okay, now I need to calculate the angular velocity vector associated
        % with each wing.  This will provide the stroke deviation. Calculate
        % these from the quaternion values
        % q_dot = .5 q * omega_F
        % q_dot = .5 omega_B * q


        %===============================
        % Right wing
        %===============================

        quat = SOLN(12:15,:);

        %quat = [-quat(1:3,:) ; quat(4,:)];

        %Calculate the wing tip trajectory & leading edge orientation
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
        % Left Wing
        quat = SOLN(8:11,:);
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


        %%=======================
        % Calculate the body orientation
        for k = 1:length(movidx)
            Rbod = quat2matNEW(SOLN(4:7,k));
            BodyAng_auto(:,k) = Rot2Euler(Rbod');%Rot2Euler(G(1:3,1:3)');
        end


        %========================
        % Apply an extra rotation about the fixed body axis to calculate the
        % wing rotations in the appropriate stroke plane
        % NomPitch =  acos(vv(1))*(180/pi);

        % Set the pitch fixed to conform with that measured in Charlie David paper for
        % hovering.
        NomPitch = 62;

        q_mod = [0; sin(0.5.*NomPitch.*(pi/180)) ; 0; cos(0.5*NomPitch.*(pi/180))];

        %Save this Rotation in the Manual Fit Structure
        ManualFit.NomPitch = NomPitch;
        ManualFit.q_mod = q_mod;
        save([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub],'ManualFit');

        %%================================
        % Add this rotation about the y axis to the wing quaternions so that 
        % the mean stroke plane aligns with the horizontal.
        Rquat = zeros(size(SOLN(12:15,:)));
        Lquat = zeros(size(SOLN(8:11,:)));

        for k = 1:size(Rquat,2)
            Rquat(:,k) = quatprod(q_mod,SOLN(12:15,k));
            Lquat(:,k) = quatprod(q_mod,SOLN(8:11,k));
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

        %================================
        % Now, smooth the wing and body angles
        % Cutoff frequency from fft on a sample signal
        filter_order = 4;
        [b a] = butter( filter_order,Fcutoff_body*(2/Fs));

        for k = 1:size(BodyAng_auto,1)
            BodyAng_auto(k,:) = filtfilt(b,a,BodyAng_auto(k,:));
        end

        [b a] = butter( filter_order,Fcutoff_wing*(2/Fs));

        phi_R = filtfilt(b,a,phi_R);
        phi_L = filtfilt(b,a,phi_L);

        theta_R = filtfilt(b,a,theta_R);
        theta_L = filtfilt(b,a,theta_L);

        alpha_R = filtfilt(b,a,alpha_R);
        alpha_L = filtfilt(b,a,alpha_L);   

    end
end