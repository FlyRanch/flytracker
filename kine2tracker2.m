clc; 
%close all; 
clear all
load('ICFromKine/cisco098.mat');
%%
load flygenmod
params = flygenmod;

startpath = '/Users/cisco/FlyVideo/solutions/';

frame = 325;
frmss = 325:425;

PAR = LoadVideo;
PAR.c = 4;
PAR.CalibFile = 'kine/Calibrations for Kine/cal_for_20031217_done_20040205';
%
%
%
%
% - Camera Info
PAR.dt = 1/6000;  %Framerate of the camera
PAR.numcam = 3;
PAR.imgres = [512 512];
%This is the number of sample points used within the Fly model.
%It changes how fine the mesh is. 
PAR.L1 = 15; %# of steps for body along length
PAR.L2 = 6; %# of steps for head along length
PAR.L3 = 30; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

PAR.flynum = [1 ];  %[1 2];
PAR.numfly = length(PAR.flynum);

%The subsample scale of the image is 1/2^(PAR.pwr)
PAR.etamax = 0;
PAR.paramdim = 1;
PAR.statedim = 15*ones(1,PAR.numfly); 
PAR.pNoisedim = PAR.statedim;

PAR.OccludeShape = {[],[],[]};

PAR.NumPrevSol = 1;

PAR.streams = 2;
        
ManualFit.ImageData = PAR;


if exist([PAR.solutionpath ['fly_' PAR.stub]],'dir') ~= 7
    mkdir(PAR.solutionpath,['fly_' PAR.stub]);
    mkdir(PAR.solutionpath,['Features_' PAR.stub]);
end

save([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub],'ManualFit');

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

for iii = 1:100

    %===============================================================
    %% Get body params and coordinates
    %params.bodyscale = data.kine.body.data.length(frmss(iii)) /
    %params.body_tip2tip;
    params.bodyscale = data.kine.body.data.length(frmss(iii)) / (params.bodylen+params.headlen);
    avg_winglength = mean([data.kine.left_wing.data.length(frmss(iii))
    data.kine.right_wing.data.length(frmss(iii))]);
    params.wingscale = avg_winglength / params.wing_tip2tip;

    %===============================================================
    %% Get the body orientation
    q_body = data.kine.body.data.quat(:,frmss(iii));

    %I store the scalar part of quaternion at the end
    q_body = [q_body(2:4); q_body(1)];

    % add a rotation by pi along the roll axis since our model has the body
    % fixed frame with z-axis pointing ventral.
    q_body = quatprod(q_body,[1*sin(pi/2) 0 0 cos(pi/2)]);

    T_body = data.kine.body.data.v_trans(:,frmss(iii));

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


    %==============================================================
    %
    %
    %
    %===============================================================
    % Get the wing orientations
    %% Left
    q_Lwing = data.kine.left_wing.data.quat(:,frmss(iii));

    %I store the scalar part of quaternion at the end
    q_Lwing = [q_Lwing(2:4); q_Lwing(1)];

    % Premultiply rotation by this alignment quaternion that takes into account
    % the orientation of coordinate axis fixed to the left wing.
    q_Lwingaxisalign = quat2matNEW([0 -1 0;1 0 0;0 0 1]);
    q_Lwing = quatprod(q_Lwing,q_Lwingaxisalign);


    T_L = data.kine.left_wing.data.v_trans(:,frmss(iii));
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
    q_Rwing = data.kine.right_wing.data.quat(:,frmss(iii));

    %I store the scalar part of quaternion at the end
    q_Rwing = [q_Rwing(2:4); q_Rwing(1)];

    % Premultiply rotation by this alignment quaternion that takes into account
    % the orientation of coordinate axis fixed to the right wing.
    q_Rwingaxisalign = quat2matNEW([0 1 0;1 0 0;0 0 -1]);
    q_Rwing = quatprod(q_Rwing,q_Rwingaxisalign);


    T_R = data.kine.right_wing.data.v_trans(:,frmss(iii));
    T_R = T_R(:);

    %Just take the orientation because I assume that the wing is fixed at the
    %joint.
    q_Rwing_rel = quatprod([-q_body(1:3); q_body(4)],q_Rwing);

    % Relative translation from body fixed axis to wing joint
    params.T_Rwing_rel = qxform([-q_body(1:3) ; q_body(4)],T_R - T_body);

    %-----------------------------------------------
    %% Assemble quaternion based fly state
    %pQ = [T_body
    %    q_body
    %    q_Lwing_rel
    %    q_Rwing_rel];

    % If the fourth element in the wing quaternions is not
    % positive, then change the sign of the the wing quaternions 
    % to make it match with the database prior.  Remember, q and -q
    % map to the same point in SO(3).
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

    ManualFit.frame = frmss(iii);
    ManualFit.soln = p;
    ManualFit.solnQ = pQ;
    ManualFit.DLT = DLT;
    ManualFit.params = params;
    ManualFit.cam = cam;

    %Plot initial condition to check for correctness

    paste_imagefunQ(pQ,frmss(iii),ManualFit,PAR);
    if iii == 1
        save([PAR.solutionpath 'fly_' PAR.stub '/ManualFit_' PAR.stub],'ManualFit');
    end
    

display('tee-hee')
pause
end