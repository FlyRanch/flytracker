function pnew = fixRollTPlaneQ(p,PAR)

%Euler angle rotation functions
Rx = @(theta) [1 0 0
    0 cos(theta) -sin(theta)
    0 sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 sin(theta)
    0 1 0
    -sin(theta) 0  cos(theta)];
Rz = @(theta) [cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];

%Convert quaternion state to the joint
Tbody = p(1:3);
qbody = p(4:7);

BodyOnly = length(p) < 8;
Velstate = length(p) > 15;
if ~BodyOnly
    qLwing = p(8:11);
    qRwing = p(12:15);
end
if Velstate
    Vel = p(16:21);
end
R_body = quat2matNEW(qbody);

Gbody = [R_body Tbody ; zeros(1,3) 1];
% S_body = homo2screw(Gbody);
% 
% Theta_Lwing = Rot2Joint(quat2matNEW(qLwing));
% Theta_Rwing = Rot2Joint(quat2matNEW(qRwing));
% 
% p = [S_body
%     Theta_Lwing
%     Theta_Rwing];
% 
% 
% % Grab spatial transformations
% [Gbody,GL,GR] = GTransform(p);

GL = quat2matNEW(quatprod(qbody,qLwing));
GR = quat2matNEW(quatprod(qbody,qRwing));

%Wing vectors that point from wing tip to wing hinge
Vr = -GR(1:3,2);
Vl = GL(1:3,2);

%Coordinate frame attached to body
% Zb = Gbody(1:3,3);
% Yb = Gbody(1:3,2);
% Xb = Gbody(1:3,1);
% Tb = Gbody(1:3,4);
Zb = R_body(1:3,3);
Yb = R_body(1:3,2);
Xb = R_body(1:3,1);
Tb = Tbody;

%wing Vectors in the transverse plane defined by [Zb Yb]
Vlbar = [Vl'*Zb Vl'*Yb];
Vrbar = [Vr'*Zb Vr'*Yb];

%if the norm of these vectors is small, the wings are tucked back, and are
% not extended outward, so the projection of
%the wing vectors (Vlbar, Vrbar) is very small.  If this is the case, do
%not perform the roll update because the fly is not flying
Wingflag = (norm(Vlbar) < 0.4 || norm(Vrbar) < 0.4);

%Normalize
Vrbar = Vrbar./norm(Vrbar);
Vlbar = Vlbar./norm(Vlbar);


%keyboard
% move Vr into first quadrant and take inner product
% the acos is the angle between the two vectors
tt = acos([Vrbar(1) -Vrbar(2)]*Vlbar');

if acos(Vrbar(1)) < acos(Vlbar(1))
    tt = -tt;
elseif acos(Vrbar(1)) > acos(Vlbar(1))
    tt = tt;
end

%keyboard
%If the deviation is too large, don't apply the roll correction because the
%animal is probably not airborne and has intentional asymmetry in the wing
%locations

if abs((tt/2)*(180/pi)) > 80
    pnew = p;

%     %Convert joint representation back to quaternion state
%     bodyscrew = pnew(1:6);
%     Gbody = screw2homo(bodyscrew);
% 
%     Tbody = Gbody(1:3,4);
%     Rbody = Gbody(1:3,1:3);
%     q_body = quat2matNEW(Rbody);
% 
%     LTheta = pnew(7:9);
%     RTheta = pnew(10:12);
% 
%     qL = quat2matNEW(Rx(LTheta(1))*Ry(LTheta(2))*Rz(LTheta(3)));
%     qR = quat2matNEW(Rx(RTheta(1))*Ry(RTheta(2))*Rz(RTheta(3)));
% 
%     pnew = [Tbody
%         q_body
%         qL
%         qR];

    fprintf(['No Roll Correction\n']);
    return
elseif Wingflag
    pnew = p;

%     %Convert joint representation back to quaternion state
%     bodyscrew = pnew(1:6);
%     Gbody = screw2homo(bodyscrew);
% 
%     Tbody = Gbody(1:3,4);
%     Rbody = Gbody(1:3,1:3);
%     q_body = quat2matNEW(Rbody);
% 
%     LTheta = pnew(7:9);
%     RTheta = pnew(10:12);
% 
%     qL = quat2matNEW(Rx(LTheta(1))*Ry(LTheta(2))*Rz(LTheta(3)));
%     qR = quat2matNEW(Rx(RTheta(1))*Ry(RTheta(2))*Rz(RTheta(3)));
% 
%     pnew = [Tbody
%         q_body
%         qL
%         qR];

    fprintf(['No Roll Correction\n']);
    return
else
    %fprintf(['Roll Correction\n']);
    %Apply this rotation about the body X-axis (roll);
    Groll = screw2homo([cross(-sign(tt)*Xb,Tb);(tt/2).*Xb]);

    %Store new body translation;
    pnew = p;
    %pnew(1:6) = homo2screw(Groll*Gbody);
    Gbodynew = Groll*Gbody;
    pnew(1:3) = Gbodynew(1:3,4);
    pnew(4:7) = quat2matNEW(Gbodynew(1:3,1:3));
    
    %Now, update the wing joint parameters so that they match the original
    %location

    clear flymod
    [x,y,z] = flymod(p,PAR.params,PAR);

    pts = cell(3,1);
    %Grab the wing points because I want keep them fixed
    for j = 2:length(x)
        tmp = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
        pts{j,:} = reshape(tmp',[],1);
    end
    WingPts = cell2mat(pts);

    pfull = pnew;
    p0 = pfull(7:12);
    options = optimset('display','none');
    JointNew = lsqnonlin(@PtErr,p0,[],[],options);

    pnew(7:12) = JointNew;

    %Convert joint representation back to quaternion state
    bodyscrew = pnew(1:6);
    Gbody = screw2homo(bodyscrew);

    Tbody = Gbody(1:3,4);
    Rbody = Gbody(1:3,1:3);
    q_body = quat2matNEW(Rbody);

    LTheta = pnew(7:9);
    RTheta = pnew(10:12);

    qL = quat2matNEW(Rx(LTheta(1))*Ry(LTheta(2))*Rz(LTheta(3)));
    qR = quat2matNEW(Rx(RTheta(1))*Ry(RTheta(2))*Rz(RTheta(3)));

    pnew = [Tbody
        q_body
        qL
        qR];
    
    if Velstate
        % update the angular velocities based on the roll correction factor
        % Instantaneous velocity based on twist
        twist = [cross(-sign(tt)*Xb,Tb) ; Xb];
        %this is the theta-dot, the amount of rotation per unit time
        theta_dot = (tt/2) / PAR.dt;
        
        %angular displacement
        disp = Vel(4:6).*PAR.dt;
        %add the displacement due to the roll correction
        disp = disp + (tt/2)*twist(4:6);
        %New velocity is this displacement per unit time
        Vel(4:6) = disp./PAR.dt;        
        
        pnew = [pnew ; Vel];
    end
    
end

    function [Gbody,GL,GR] = GTransform(p)

        BL = PAR.params.bodyscale*(PAR.params.bodylen+PAR.params.headlen);
        RJTrans = BL.*([0.2021 0.1055 -0.1477] );
        LJTrans = BL.*([0.2021 -0.1055 -0.1477]);


        bodyscrew = p(1:6);
        Gbody = screw2homo(bodyscrew);
        Tbody = Gbody(1:3,4);
        Rbody = Gbody(1:3,1:3);



        Theta = p(7:12);
        S = [zeros(3,6) ; [eye(3) eye(3)]];
        %Scale the twists according to joint angle
        S(4:6,:) = S(4:6,:).*repmat(Theta',3,1);

        G_BL = screw2homo(S(:,1))*screw2homo(S(:,2))*screw2homo(S(:,3));
        G_BR = screw2homo(S(:,4))*screw2homo(S(:,5))*screw2homo(S(:,6));


        % Calculate full transformation
        G_BL(1:3,1:3) = G_BL(1:3,1:3)*PAR.params.wingscale;
        G_BL(1:3,4) = LJTrans';

        G_BR(1:3,1:3) = G_BR(1:3,1:3)*PAR.params.wingscale;
        G_BR(1:3,4) = RJTrans';

        GL = Gbody*G_BL;
        GR = Gbody*G_BR;

    end

    function F = PtErr(pp)
        persistent xorg yorg zorg



        x = cell(2,1);
        y = x;
        z = x;
        pts = x;

        if isempty(xorg)
            %If these points don't exist from a previous function call, create them

            %=========================================================
            %WING
            %=========================================================
            %[x{2},y{2},z{2}] = flywingmod(PAR.params.wing,PAR.params.winglen,PAR);
            [x{1},y{1},z{1}] = flywingmodL(PAR.params.wing,PAR.params.winglen,PAR);

            %Make Right Wing too
            % x{3} = x{2};
            % y{3} = y{2};
            % z{3} = z{2};
            [x{2},y{2},z{2}] = flywingmodR(PAR.params.wing,PAR.params.winglen,PAR);

            xorg = x;
            yorg = y;
            zorg = z;
        else
            % these points already exist, so just assign them to the appropriate
            % variables
            x = xorg;
            y = yorg;
            z = zorg;
        end

        BL = PAR.params.bodyscale*(PAR.params.bodylen+PAR.params.headlen);
        RJTrans = BL.*([0.2021 0.1055 -0.1477] );
        LJTrans = BL.*([0.2021 -0.1055 -0.1477]);

        %Assemble full state

        p = pfull;
        p(7:12) = pp;

        bodyscrew = p(1:6);
        Gbody = screw2homo(bodyscrew);
        Tbody = Gbody(1:3,4);
        Rbody = Gbody(1:3,1:3);


        Theta = p(7:12);
        S = [zeros(3,6) ; [eye(3) eye(3)]];
        %Scale the twists according to joint angle
        S(4:6,:) = S(4:6,:).*repmat(Theta',3,1);

        G_BL = screw2homo(S(:,1))*screw2homo(S(:,2))*screw2homo(S(:,3));
        G_BR = screw2homo(S(:,4))*screw2homo(S(:,5))*screw2homo(S(:,6));

        %The translation part of the homogeneous transformations above will always
        %be zero, so just include the wing joint translations below
        %keyboard

        [x{1},y{1},z{1}] = xform_surf(x{1},y{1},z{1},LJTrans,G_BL(1:3,1:3),PAR.params.wingscale);
        [x{2},y{2},z{2}] = xform_surf(x{2},y{2},z{2},RJTrans,G_BR(1:3,1:3),PAR.params.wingscale);

        for k = 1:length(x)
            [x{k},y{k},z{k}] = xform_surf(x{k},y{k},z{k},Tbody,Rbody);
        end

        for i = 1:length(x)
            tmp = [reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
            pts{i,:} = reshape(tmp',[],1);
        end
        NewWingPts = cell2mat(pts);


        F = NewWingPts - WingPts;

    end

end