% FOR is right handed, with -z used only for display, angles are anticlockwise positive (looking along the axis down onto origin)
% Calculations are performed for stroke plane

function calcBodyStrokePlaneCentred(theFrame)

if theFrame == 259
    'hello';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Angles in common frame of reference (right handed, anticlockwise +), in radians
% yaw, pitch and roll are in a body centred FOR (as opposed to stroke plane FOR)
B_yaw = 0;         % AC+, looking down onto fly, zero facing x/0/0, 0-360 deg
B_pitch = 0;       % Positive down, zero horizontal, 0-360 deg
B_roll = 0;        % C+ in heading direction, zero horizontal (this is because model fly looks along x axis, and AC+ applies)

% Right wing
% Wing positions are relative to the stroke plane!
rightAzi = 0;       % AC+ relative to x axis
rightEle = 30;       % anticlockwise +, i.e. positive angles depress the wing
rightAlpha = 0;     % 0 is stroke plane, anticlockwise positive looking down from wing tip to hinge!!!

% Left wing
leftAzi = 0;        % AC+ relative to x axis
leftEle = 0;        % anticlockwise +, i.e. positive angles depress the wing
leftAlpha = 0;      % 0 is stroke plane, anticlockwise positive looking down from wing tip to hinge!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Button values are:
head=1;tail=2;rwh=3;rwt=4;lwh=5;lwt=6;rwn=7;lwn=8;
right=1;left=2; 
x=1;y=2;z=3;

global KineData;

%theFrame = KineData.FrameSliderValue;

% I CHANGE SIGN OF Z NOW. THIS PUTS FLY INTO A X,Y,Z FRAME OF REFERENCE. -Z IS USED ONLY FOR DISPLAY. 
% Put all body positions into 1 vector. Change sign on z now, to use right handed coord frame x,y,z. 
v_head=[KineData.pos.x(head,theFrame);KineData.pos.y(head,theFrame);-KineData.pos.z(head,theFrame)];
v_tail=[KineData.pos.x(tail,theFrame);KineData.pos.y(tail,theFrame);-KineData.pos.z(tail,theFrame)];
v_rwh=[KineData.pos.x(rwh,theFrame);KineData.pos.y(rwh,theFrame);-KineData.pos.z(rwh,theFrame)];
v_rwt=[KineData.pos.x(rwt,theFrame);KineData.pos.y(rwt,theFrame);-KineData.pos.z(rwt,theFrame)];
v_lwh=[KineData.pos.x(lwh,theFrame);KineData.pos.y(lwh,theFrame);-KineData.pos.z(lwh,theFrame)];
v_lwt=[KineData.pos.x(lwt,theFrame);KineData.pos.y(lwt,theFrame);-KineData.pos.z(lwt,theFrame)];

% Calc right wing normal
dx_r=v_rwt(1) - v_rwh(1);                                 % The next 6 lines of code come directly from updateWingModel
dy_r=v_rwt(2) - v_rwh(2);
dz_r=v_rwt(3) - v_rwh(3);
[theta_r,phi_r,length_r] = cart2sph(dx_r,dy_r,dz_r); % anticlockwise + in right handed coord frame

alpha_r = KineData.pos.angles(right,theFrame)-pi;
z_dash_dash_r=rot3D([0; 0; 1],-theta_r,-phi_r,-alpha_r); % direction of vector is up when looking from the tip down onto the hinge with the leading edge on the right
v_rwn = z_dash_dash_r + v_rwh;
KineData.pos.x(rwn,theFrame)=z_dash_dash_r(1);
KineData.pos.y(rwn,theFrame)=z_dash_dash_r(2);
KineData.pos.z(rwn,theFrame)=z_dash_dash_r(3);


% Calc left wing normal
dx_l=v_lwt(1) - v_lwh(1);
dy_l=v_lwt(2) - v_lwh(2);
dz_l=v_lwt(3) - v_lwh(3);
[theta_l,phi_l,length_l] = cart2sph(dx_l,dy_l,dz_l);

alpha_l = -(KineData.pos.angles(left,theFrame)-pi);
z_dash_dash_l=rot3D([0; 0; 1],-theta_l,-phi_l,alpha_l);% direction of vector is up when looking from the tip down onto the hinge with the leading edge on the right
v_lwn = z_dash_dash_l + v_lwh;
KineData.pos.x(lwn,theFrame)=z_dash_dash_l(1);
KineData.pos.y(lwn,theFrame)=z_dash_dash_l(2);
KineData.pos.z(lwn,theFrame)=z_dash_dash_l(3);% wing normal in right handed frame of reference!!! 


v_fly=[v_head v_tail v_rwh v_rwt v_lwh v_lwt v_rwn v_lwn];

% Shift fly so that tail is at 0,0,0
v_fly(1,:)=v_fly(1,:)-v_tail(1,:); % x. normalise x of all points to tail position = 0
v_fly(2,:)=v_fly(2,:)-v_tail(2,:); % y
v_fly(3,:)=v_fly(3,:)-v_tail(3,:); % z

% Rotate the fly such that the head lies on the x axis
dx = v_head(1) - v_tail(1);
dy = v_head(2) - v_tail(2);
dz = v_head(3) - v_tail(3);
[theta,phi,length] = cart2sph(dx,dy,dz);                                            % Position of head relative to the tail in Spherical coords. 

% stroke plane axis is 60 deg different in phi
%phi = -(phi - pi/3); % 60 deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_pitch = -phi;       % Positive down, zero horizontal, 0-360 deg
B_yaw = theta;         % AC+, looking down onto fly, zero facing x/0/0, 0-360 deg
%test_SP_yaw = SP_yaw / pi*180
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_fly_bodyCentred=invRot3D(v_fly,-theta,-phi,0);                                    % [THETA,PHI,R] = cart2sph(X,Y,Z); invRot3D(v_fly,-THETA,-PHI,0); This always yields x,0,0. x = len(X,Y,Z)



% This algo to roll the fly to the horizontal using the 3D wing hinge positions
% Rolling has to be done after the conversion into a body centred frame of reference.
%aRoll=atan2(v_fly_bodyCentred(3,3)-v_fly_bodyCentred(3,5),v_fly_bodyCentred(2,3)-v_fly_bodyCentred(2,5));       % zRWH-zLWH, yRWH-zLWH

dx_hinge_BC=v_fly_bodyCentred(1,rwh)-v_fly_bodyCentred(1,lwh);
dy_hinge_BC=v_fly_bodyCentred(2,rwh)-v_fly_bodyCentred(2,lwh);
dz_hinge_BC=v_fly_bodyCentred(3,rwh)-v_fly_bodyCentred(3,lwh);
[THETA,aRoll,R] = cart2sph(dy_hinge_BC,dx_hinge_BC,dz_hinge_BC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_roll = -aRoll;        % C+ in heading direction, zero horizontal (this is because model fly looks along x axis, and AC+ applies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_fly_bodyCentred=invRot3D(v_fly_bodyCentred,0,0,-aRoll);

% Calc body centred for display only
calcBodyCentred(v_fly_bodyCentred,theFrame);

% STROKE PLANE CENTRED BEGINS HERE. We start with v_fly, with the tail in 0,0,0
v_fly_strokePlaneCentred=invRot3D(v_fly_bodyCentred,0,KineData.StrokePlaneDeg/180*pi,0);    % Rotation around Y axis, i.e. tail

% RIGHT WING
% Calc wing angle relative to horizontal plane, i.e. body plane

rwn_relative = v_fly_strokePlaneCentred(:,rwn) - v_fly_strokePlaneCentred(:,rwh);        % full cicle
r_wing_tip_relative = v_fly_strokePlaneCentred(:,rwt) - v_fly_strokePlaneCentred(:,rwh);
dx_r = r_wing_tip_relative(1);
dy_r = r_wing_tip_relative(2);
dz_r = r_wing_tip_relative(3);

% The following 2 steps would rotate the wing tip onto the y axis. I rotate rwn_relative instead, corresponding to the wing at 0 ele and azi. 
% 1. rotate wing tip such that it lies on YZ Plane (x=0)
theta_r = atan2(dy_r,dx_r);% angle between x and [dx dy]
%theta = theta +pi/2;% angle between y and [dx dy] % AC+, 0 is (0/-1/0)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rightAzi=theta_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotate wing tip to point along x axis, apply to wing normal
rwnXZPlane=arbitraryRotate(rwn_relative,-theta_r,[0 0 1]);% YZ refers to tip postion
test1_r=arbitraryRotate([dx_r;dy_r;dz_r],-theta_r,[0 0 1]);



% 2. rotate wing tip such that it lies on x axis
phi_r = atan2(dz_r,sqrt(dx_r^2+dy_r^2));% elevation
rwnXAxis=arbitraryRotate(rwnXZPlane,phi_r,[0 1 0]);%X refers to tip postion
test2_r=arbitraryRotate(test1_r,phi_r,[0 1 0]);%X refers to tip postion
%phi = -phi; % AC+ makes elevation negative, depression positive. That's a reason for hope...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rightEle=phi_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wingAngle_r = atan2(rwnXAxis(3),rwnXAxis(2));
wingAngle_r = +(wingAngle_r-pi/2);                                              %anticlockwise + when looking from tip to hinge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rightAlpha=wingAngle_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test_rwn_relative = rot3D([0; 0; 1],-pi/2+theta,-phi,wingAngle_r);
%test_rwn_relative(1)= -test_rwn_relative(1)
%rwn_relative - test_rwn_relative;

% LEFT WING
lwn_relative = v_fly_strokePlaneCentred(1:3,lwn) - v_fly_strokePlaneCentred(1:3,lwh);
l_wing_tip_relative = v_fly_strokePlaneCentred(:,lwt) - v_fly_strokePlaneCentred(:,lwh);
dx_l = l_wing_tip_relative(1);
dy_l = l_wing_tip_relative(2);
dz_l = l_wing_tip_relative(3);

theta_l = atan2(dy_l,dx_l);% angle between x and [dx dy]
%theta = theta -pi/2;% angle between y and [dx dy]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leftAzi=theta_l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lwnXZPlane=arbitraryRotate(lwn_relative,-theta_l,[0 0 1]);% YZ refers to tip postion
test1_l=arbitraryRotate([dx_r;dy_r;dz_r],-theta_l,[0 0 1]);

% 2. rotate wing tip such that it lies on y axis (x=0,z=0)
phi_l = atan2(dz_l,sqrt(dx_l^2+dy_l^2));% elevation
lwnXAxis=arbitraryRotate(lwnXZPlane,phi_l,[0 1 0]);%Y refers to tip postion
test2_l=arbitraryRotate(test1_l,phi_l,[0 1 0]);%X refers to tip postion
% phi = phi; % AC+ makes elevation +

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leftEle=phi_l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wingAngle_l = atan2(lwnXAxis(3),lwnXAxis(2));
wingAngle_l = +(wingAngle_l-pi/2);                                              %anticlockwise + when looking from tip to hinge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leftAlpha=wingAngle_l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the data into the global struct
KineData.posBCStrokePlane.x(head,theFrame) = v_fly_strokePlaneCentred(1,head);
KineData.posBCStrokePlane.y(head,theFrame) = v_fly_strokePlaneCentred(2,head);
KineData.posBCStrokePlane.z(head,theFrame) = -v_fly_strokePlaneCentred(3,head);

KineData.posBCStrokePlane.x(tail,theFrame) = v_fly_strokePlaneCentred(1,tail);
KineData.posBCStrokePlane.y(tail,theFrame) = v_fly_strokePlaneCentred(2,tail);
KineData.posBCStrokePlane.z(tail,theFrame) = -v_fly_strokePlaneCentred(3,tail);

KineData.posBCStrokePlane.x(rwh,theFrame) = v_fly_strokePlaneCentred(1,rwh);
KineData.posBCStrokePlane.y(rwh,theFrame) = v_fly_strokePlaneCentred(2,rwh);
KineData.posBCStrokePlane.z(rwh,theFrame) = -v_fly_strokePlaneCentred(3,rwh);

KineData.posBCStrokePlane.x(lwh,theFrame) = v_fly_strokePlaneCentred(1,lwh);
KineData.posBCStrokePlane.y(lwh,theFrame) = v_fly_strokePlaneCentred(2,lwh);
KineData.posBCStrokePlane.z(lwh,theFrame) = -v_fly_strokePlaneCentred(3,lwh);

KineData.posBCStrokePlane.x(rwt,theFrame) = v_fly_strokePlaneCentred(1,rwt);
KineData.posBCStrokePlane.y(rwt,theFrame) = v_fly_strokePlaneCentred(2,rwt);
KineData.posBCStrokePlane.z(rwt,theFrame) = -v_fly_strokePlaneCentred(3,rwt);

KineData.posBCStrokePlane.x(lwt,theFrame) = v_fly_strokePlaneCentred(1,lwt);
KineData.posBCStrokePlane.y(lwt,theFrame) = v_fly_strokePlaneCentred(2,lwt);
KineData.posBCStrokePlane.z(lwt,theFrame) = -v_fly_strokePlaneCentred(3,lwt);

% HERE THE WING ANGLE REL TO THE STROKE PLANE IS CALCULATED
KineData.posBCStrokePlane.angles(left,theFrame) = wingAngle_l;
KineData.posBCStrokePlane.angles(right,theFrame) = wingAngle_r;


% Robofly angles in degrees. Wing coordinates are relative to stroke plane (defined as 60 deg rotated around pitch axis). 

% Left wing
roboAzi_l = +(leftAzi/pi*180 - 90);        % frontal positions positive, lateral = 0
roboEle_l = leftEle/pi*180;        % 0 is horzontal, + up
roboAlpha_l  = +(leftAlpha/pi*180 + 90) ;    % Positive angle of attack on downstroke

% Right wing
roboAzi_r = -(rightAzi/pi*180 + 90);;        % frontal positions positive, lateral = 0
roboEle_r = rightEle/pi*180;        % 0 is horzontal, + up
roboAlpha_r  = -(rightAlpha/pi*180 + 90);    % Positive angle of attack on downstroke

yaw = B_yaw/pi*180;       % left turn +, i.e. AC+, looking down onto fly, zero facing x/0/0, 0-360 deg
pitch = B_pitch/pi*180;      % BODY pitch (not SP). We have already rotated the fly by 60 deg, so subtract that. Positive down, zero horizontal, 0-360 deg
roll = B_roll/pi*180;       % C+ in heading direction, zero horizontal (this is because model fly looks along x axis, and AC+ applies)

%%% Calculate Rotation Matrix for ZYX rotation from above values

% Note that the X axis points from the tail to the head!!!! 

t = B_yaw;          % left turn +, i.e. AC+, looking down onto fly
ph = B_pitch;      % Positive DOWN, i.e. AC+ looking from left wing to right wing
ps = B_roll;       % AC+ looking from head to tail

rx=[[1,0,0];[0, cos(ps), -sin(ps)];[0,sin(ps),cos(ps)]];
ry=[[cos(ph),0,sin(ph)];[0,1,0];[-sin(ph),0,cos(ph)]];
rz=[[cos(t),-sin(t),0];[sin(t),cos(t),0];[0,0,1]];

if theFrame == 100
    'hello';
end

R = rz * ry * rx; % This is the rotation matrix according to Haslwanter, and also http://www.magic-software.com/Documentation/EulerAngles.pdf
KineData.bodyMotion.R(1:3,1:3,theFrame) = R;% Rotation matrix using ZYX Euler Angles
% Factoring rotation matrix R according with sequence YXZ. thetaZ is morphological yaw!  
r20=R(3,1);r21=R(3,2);r22=R(3,3);
r10=R(2,1);r11=R(2,2);r12=R(2,3);
r00=R(1,1);r01=R(1,2);r02=R(1,3);

thetaX = asin(-r12);
if thetaX < pi/2
    if thetaX > -pi/2
        thetaY = atan2(r02,r22);
        thetaZ = atan2(r10,r11);
    else
        % not a unique solution
        thetaY = -atan2(-r01,r00);
        thetaZ = 0;
    end
else    % not a unique solution
    thetaY = atan2(-r01,r00);
    thetaZ = 0;
end

%thetaX/pi*180;
%thetaY/pi*180;
thetaZ = thetaZ/pi*180;




KineData.robo.theta.left(theFrame) =  angleAroundZeroDegrees(roboEle_l);
KineData.robo.theta.right(theFrame) =  angleAroundZeroDegrees(roboEle_r);


KineData.robo.phi.left(theFrame) =  angleAroundZeroDegrees(roboAzi_l);
KineData.robo.phi.right(theFrame) =  angleAroundZeroDegrees(roboAzi_r);

KineData.robo.alpha.left(theFrame) =  angleAroundZeroDegrees(roboAlpha_l);
KineData.robo.alpha.right(theFrame) =  angleAroundZeroDegrees(roboAlpha_r);

KineData.bodyMotion.yaw(theFrame) = thetaZ;   % AC+ looking down onto fly
KineData.bodyMotion.pitch(theFrame) = pitch; % positive down, i.e. clockwise - if looking down the left wing
KineData.bodyMotion.roll(theFrame) = roll;   % right bank +, i.e. AC+ looking from front


% For the MFI people, I also need body coordinates in XYZ sequence....

thetaY = asin(r02);
if thetaY < pi/2
    if thetaY > -pi/2
        thetaX = atan2(-r12,r22);
        thetaZ = atan2(-r01,r00);
    else
        % not a unique solution
        thetaX = -atan2(r10,r11);
        thetaZ = 0;
    end
else
    % not a unique solution
    thetaX = atan2(r10,r11);
    thetaZ = 0;
end

thetaX=thetaX/pi*180;
thetaY=thetaY/pi*180;
thetaZ = thetaZ/pi*180;

KineData.bodyMotion.XYZ.yaw(theFrame) = thetaZ;
KineData.bodyMotion.XYZ.pitch(theFrame) = thetaY;
KineData.bodyMotion.XYZ.roll(theFrame) = thetaX;





%************************************************************************************************************%
%************************************************************************************************************%
%************************************************************************************************************%
%************************************************************************************************************%
%************************************************************************************************************%
%************************************************************************************************************%
%************************************************************************************************************%




function calcBodyCentred(v_fly_bodyCentred,theFrame)

head=1;tail=2;rwh=3;rwt=4;lwh=5;lwt=6;rwn=7;lwn=8;
right=1;left=2;

global KineData


% RIGHT WING
% Calc wing angle relative to horizontal plane, i.e. body plane
rwn_relative = v_fly_bodyCentred(:,rwn) - v_fly_bodyCentred(:,rwh);
r_wing_tip_relative = v_fly_bodyCentred(:,rwt) - v_fly_bodyCentred(:,rwh);
dx_r = r_wing_tip_relative(1);
dy_r = r_wing_tip_relative(2);
dz_r = r_wing_tip_relative(3);

% The following 2 steps would rotate the wing tip onto the y axis. I rotate rwn_relative instead, corresponding to the wing at 0 ele and azi. 
% 1. rotate wing tip such that it lies on YZ Plane (x=0)
theta_r = atan2(dy_r,dx_r);% angle between x and [dx dy]
rwnXZPlane=arbitraryRotate(rwn_relative,-theta_r,[0 0 1]);% YZ refers to tip postion

% 2. rotate wing tip such that it lies on y axis (x=0,z=0)
phi_r = atan2(dz_r,sqrt(dx_r^2+dy_r^2));% elevation
rwnXAxis=arbitraryRotate(rwnXZPlane,phi_r,[0 1 0]);%X refers to tip postion
wingAngle_r = atan2(rwnXAxis(3),rwnXAxis(2));
wingAngle_r = +(wingAngle_r-pi/2);                                              %anticlockwise + when looking from tip to hinge


% LEFT WING
lwn_relative = v_fly_bodyCentred(:,lwn) - v_fly_bodyCentred(:,lwh);
l_wing_tip_relative = v_fly_bodyCentred(:,lwt) - v_fly_bodyCentred(:,lwh);
dx_l = l_wing_tip_relative(1);
dy_l = l_wing_tip_relative(2);
dz_l = l_wing_tip_relative(3);

theta_l = atan2(dy_l,dx_l);% angle between x and [dx dy]
lwnXZPlane=arbitraryRotate(lwn_relative,-theta_l,[0 0 1]);% YZ refers to tip postion

% 2. rotate wing tip such that it lies on y axis (x=0,z=0)
phi_l = atan2(dz_l,sqrt(dx_l^2+dy_l^2));% elevation
lwnXAxis=arbitraryRotate(lwnXZPlane,phi_l,[0 1 0]);%Y refers to tip postion
wingAngle_l = atan2(lwnXAxis(3),lwnXAxis(2));
wingAngle_l = +(wingAngle_l-pi/2);                                              %anticlockwise + when looking from tip to hinge




KineData.posBC.angles(left,theFrame) = wingAngle_l;
KineData.posBC.angles(right,theFrame) = wingAngle_r;

% Read the data into the global struct
KineData.posBC.x(head,theFrame) = v_fly_bodyCentred(1,head);
KineData.posBC.y(head,theFrame) = v_fly_bodyCentred(2,head);
KineData.posBC.z(head,theFrame) = -v_fly_bodyCentred(3,head);

KineData.posBC.x(tail,theFrame) = v_fly_bodyCentred(1,tail);
KineData.posBC.y(tail,theFrame) = v_fly_bodyCentred(2,tail);
KineData.posBC.z(tail,theFrame) = -v_fly_bodyCentred(3,tail);

KineData.posBC.x(rwh,theFrame) = v_fly_bodyCentred(1,rwh);
KineData.posBC.y(rwh,theFrame) = v_fly_bodyCentred(2,rwh);
KineData.posBC.z(rwh,theFrame) = -v_fly_bodyCentred(3,rwh);

KineData.posBC.x(lwh,theFrame) = v_fly_bodyCentred(1,lwh);
KineData.posBC.y(lwh,theFrame) = v_fly_bodyCentred(2,lwh);
KineData.posBC.z(lwh,theFrame) = -v_fly_bodyCentred(3,lwh);

KineData.posBC.x(rwt,theFrame) = v_fly_bodyCentred(1,rwt);
KineData.posBC.y(rwt,theFrame) = v_fly_bodyCentred(2,rwt);
KineData.posBC.z(rwt,theFrame) = -v_fly_bodyCentred(3,rwt);

KineData.posBC.x(lwt,theFrame) = v_fly_bodyCentred(1,lwt);
KineData.posBC.y(lwt,theFrame) = v_fly_bodyCentred(2,lwt);
KineData.posBC.z(lwt,theFrame) = -v_fly_bodyCentred(3,lwt);

% end of function calcBodyCentred




