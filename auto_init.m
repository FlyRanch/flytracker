function ManualFit = auto_init(ManualFit,frame,mode)

% ----------------------------------
% There are 5 modes to this function:
% 'autoradius', 'IC', and 'BG'

% 'autoradius' calculates the model shape and radius function
% from the auto generated outline.

% 'radius' calculates the model shape and radius function
% by manually clicking points

% 'calib' will calculate the pixpermm of the video sequence.

% 'BG' will calculate the background image of the video sequence.

global PAR;

%==============================
% SELECT MODE
%==============================


for j = 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch mode
        case 'BG'
            %%  Calculate segmented images
            if ~exist([PAR.solutionpath '/BG.mat'],'file')
                BG = GetBackgroundv1(PAR.ICframe + 50,PAR);
                save([PAR.solutionpath '/BG.mat'],'BG')
            else
                load([PAR.solutionpath '/BG.mat'])
            end            
            
            %Store the grayscale image backgrounds
            ManualFit.IMBG = BG;        
        case 'IC'
            load(PAR.CalibFile);
            L(:,:,1) = L_1;
            L(:,:,2) = L_2;
            L(:,:,3) = L_3;
            
            fprintf(['Refining DLT Calibration...\n']);
            for i = 1:3
                DLT(:,i) = dltfu_iter(F,L(:,:,i));
                %a(:,i) = dltfu_nonlin(F,L(:,:,i));
                cam(i) = dlt2cam(DLT(:,i));
            end

            load(['ICFromKine/' PAR.FileFromKine]);
            load flygenmod
            params = flygenmod;

            %===============================================================
            %% Get body params and coordinates
            params.bodyscale = data.kine.body.data.length(frame) / (params.bodylen+params.headlen);
            avg_winglength = mean([data.kine.left_wing.data.length(frame)
                                   data.kine.right_wing.data.length(frame)]);
            params.wingscale = avg_winglength / params.wing_tip2tip;

            %===============================================================
            %% Get the body orientation
            q_body = data.kine.body.data.quat(:,frame);

            %I store the scalar part of quaternion at the end
            q_body = [q_body(2:4); q_body(1)];

            % add a rotation by pi along the roll axis since our model has the body
            % fixed frame with z-axis pointing ventral.
            q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);

            T_body = data.kine.body.data.v_trans(:,frame);

            %Calculation translation vector -q_body*T10*q_body` + T_body
            %T10 is the approximate location of tailpt in our generative model
            %It is the next to last cross-section of model at the dorsal edge
            [xbody,ybody,zbody,s,th,X,Frenet,T10] = flybodymod(params.bodyctr,params.bodyrad,params.bodylen,PAR);

            %Say T10 is right on the tip of the thorax
            T10 = [xbody(1,1) ybody(1,1) zbody(1,1)];
            T10 = -T10;
            [T_body(1,:),T_body(2,:),T_body(3,:)] = xformq_surf(T10(1),T10(2),T10(3),T_body,q_body,1);

            %-----------------------------------------------
            %% Calculate the twist representation for these quaternions
            %             R_body = quat2matNEW(q_body);
            %
            %             G_body = [R_body T_body ; zeros(1,3) 1];
            %             S_body = homo2screw(G_body);

            %now update the body twist before calculating the relative
            %transformation from body to wings;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% First, I optimize over the body transformation so that I
            %% have the optimal localization for the fly body
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            PAR.statedim = 7;
            PAR.pNoisedim = PAR.statedim;
            PAR.BG = ManualFit.IMBG;
            PAR.params = params;
            PAR.DLT = DLT;
            PAR.cam = cam;

            %--- Initialise GSSM model from external system description script.
            model = gssm_flyOcc('init');

            % Define start and end frames of calculation
            frames = [frame frame];

            %--- Setup argument data structure which serves as input to
            %--- the 'geninfds' function. This function generates the InferenceDS
            %--- data structures which are needed by all inference algorithms
            %--- in the Rebel toolkit.

            Arg.type = 'state';                                  % inference type (state estimation)
            Arg.tag = ['State estimation for ' PAR.stub ' data.'];  % arbitrary ID tag
            Arg.model = model;                                   % GSSM data structure of external system

            % Create inference data structure and
            InfDS = geninfds(Arg);

            % generate process and observation noise sources
            [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, 'srcdkf');

            %Initialize occlusion index
            InfDS.model.Occ = cell(1,length(PAR.numfly));

            %--- initial estimate of state E[X(0)]
            %Xh(:,1) = S_body;
            Xh(:,1) = [T_body ; q_body];
            % initial state covariance
            Px_ = 0.001.*ones(1,length(Xh(:,1)));

            %Create a diagonal covariance matrix by replicating Px_ # of fish
            %times and then placing it on the diagonal of a matrix.
            Px = diag(repmat(Px_,1,PAR.numfly));

            %--- Call inference algorithm / estimator
            % Square Root Central Difference Kalman Filter
            %---------------
            InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)
            Sx = chol(Px)';

            %Update the body transformation
            [New_body, Sx, pNoise, oNoise, InternalVariablesDS,Features] = srcdkf_const(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames,1);
            %keyboard
            % Delete the initial condition file (e.g. 'fly0.mat' if frame = 1) that
            % 'srcdkf_const' saved
            
            delete([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(frame-PAR.framesample) '.mat'])

            T_body = New_body(1:3);
            q_body = New_body(4:7);
            %             % Convert back into transformation vector and quaternion
            %             Gbody = screw2homo(Sbody);
            %             T_body = Gbody(1:3,4);
            %             q_body = quat2matNEW(Gbody(1:3,1:3));

            %==============================================================
            %
            %
            %
            %===============================================================
            % Get the wing orientations
            %% Left
            q_Lwing = data.kine.left_wing.data.quat(:,frame);

            %I store the scalar part of quaternion at the end
            q_Lwing = [q_Lwing(2:4); q_Lwing(1)];

            % Premultiply rotation by this alignment quaternion that takes into account
            % the orientation of coordinate axis fixed to the left wing.
            q_Lwingaxisalign = quat2matNEW([0 -1 0;1 0 0;0 0 1]);
            q_Lwing = quatprod(q_Lwing,q_Lwingaxisalign);


            T_L = data.kine.left_wing.data.v_trans(:,frame);
            T_L = T_L(:);


            % I will calculate the relative rotations
            % for the wings by multiplying the quaternions directly.
            % Just take the orientation because I assume that the wing is fixed at the
            % joint.
            q_Lwing_rel = quatprod([-q_body(1:3) ; q_body(4)],q_Lwing);

            % Relative translation from body fixed axis to wing joint
            params.T_Lwing_rel = qxform([-q_body(1:3) ; q_body(4)],T_L - T_body);

            %--------------------------------------------
            %% Right Wing
            q_Rwing = data.kine.right_wing.data.quat(:,frame);

            %I store the scalar part of quaternion at the end
            q_Rwing = [q_Rwing(2:4); q_Rwing(1)];

            % Premultiply rotation by this alignment quaternion that takes into account
            % the orientation of coordinate axis fixed to the right wing.
            q_Rwingaxisalign = quat2matNEW([0 1 0;1 0 0;0 0 -1]);
            q_Rwing = quatprod(q_Rwing,q_Rwingaxisalign);


            T_R = data.kine.right_wing.data.v_trans(:,frame);
            T_R = T_R(:);

            %Just take the orientation because I assume that the wing is fixed at the
            %joint.
            q_Rwing_rel = quatprod([-q_body(1:3); q_body(4)],q_Rwing);

            % Relative translation from body fixed axis to wing joint
            params.T_Rwing_rel = qxform([-q_body(1:3) ; q_body(4)],T_R - T_body);

            if q_Lwing_rel(4) < 0
                q_Lwing_rel = -q_Lwing_rel;
            end
            
            if q_Rwing_rel(4) < 0
                q_Rwing_rel = -q_Rwing_rel;
            end


            pQ = [T_body
                  q_body
                  q_Lwing_rel
                  q_Rwing_rel];



            R_body = quat2matNEW(q_body);
            G_body = [R_body T_body ; zeros(1,3) 1];
            S_body = homo2screw(G_body);

            % My ordering for the twist axes corresponds to the Rxyz Euler angles.
            %
            Theta_Lwing = Rot2Joint(quat2matNEW(q_Lwing_rel));
            Theta_Rwing = Rot2Joint(quat2matNEW(q_Rwing_rel));

            p = [S_body
                Theta_Lwing
                Theta_Rwing];

            ManualFit.frame = frame;
            ManualFit.soln = p;
            ManualFit.solnQ = pQ;
            ManualFit.DLT = DLT;
            ManualFit.params = params;
            ManualFit.cam = cam;

            %Plot initial condition to check for correctness
            paste_imagefunQ(pQ,frame,ManualFit,PAR,Features);
            pause
            close all
            pause(1);
            %==============================================================         

        case 'autoradius'            
            PAR.statedim = 35 - 4;
            PAR.pNoisedim = PAR.statedim;

            %--- Initialise GSSM model from external system description script.
            model = gssm_flyShape('init');

            % Define start and end frames of calculation
            frames = [frame frame];

            Arg.type = 'state';                                  % inference type (state estimation)
            Arg.tag = ['State estimation for ' PAR.stub ' data.'];  % arbitrary ID tag
            Arg.model = model;                                   % GSSM data structure of external system

            % Create inference data structure and
            InfDS = geninfds(Arg);

            % generate process and observation noise sources
            [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, 'srcdkf');

            %Initialize occlusion index
            InfDS.model.Occ = cell(1,length(PAR.numfly));

            %--- initial estimate of state E[X(0)]
            PAR.p0 = reshape(ManualFit(1).solnQ',[],1);
            Xh = [];
            Xh(:,1) = [ManualFit.params.bodyrad(2:end-1) ; ManualFit.params.headrad(2:end-1)];
            % initial state covariance
            Px_ = 0.001.*ones(1,length(Xh(:,1)));

            %Create a diagonal covariance matrix by replicating Px_ # of fish
            %times and then placing it on the diagonal of a matrix.
            Px = diag(repmat(Px_,1,PAR.numfly));

            %--- Call inference algorithm / estimator
            % Square Root Central Difference Kalman Filter
            %---------------
            InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)
            Sx = chol(Px)';
            %keyboard
            [X,InternalVariablesDS,Features]  = srcdkf_Shape(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames);

            %keyboard
            ManualFit.params.bodyrad_old = ManualFit.params.bodyrad;
            ManualFit.params.bodyrad = [0 ; X(1:18) ; 0];

            ManualFit.params.headrad_old = ManualFit.params.headrad;
            ManualFit.params.headrad = [0 ; X(19:end) ; 0];

            %Plot initial condition to check for correctness
            paste_imagefunQ(PAR.p0,frame,ManualFit,PAR);
            pause
            close all
            pause(1);
    end
end

close all
pause(1);
