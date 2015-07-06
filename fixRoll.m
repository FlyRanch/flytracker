function p = fixRoll(p,PAR)

clear flymod
[x,y,z] = flymod(p,PAR.params,PAR);

pts = cell(3,1);
%Grab the wing points because I want keep them fixed
for j = 2:length(x)
    tmp = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
    pts{j,:} = reshape(tmp',[],1);
end
WingPts = cell2mat(pts);

Rx = @(theta) [1 0 0
    0 cos(theta) -sin(theta)
    0 sin(theta) cos(theta)];

idx = [7 9 10 12];
pfull = p;

p0 = [0 ; pfull(idx)];


%pp = lsqnonlin(@VecErr,p0);
options = optimset('display','iter','Tolfun',1e-9);
pp = fmincon(@VecErr,p0,[],[],[],[],[],[],@cfun,options);
keyboard


    function F = VecErr(pp)
        
        BL = PAR.params.bodyscale*(PAR.params.bodylen+PAR.params.headlen);
        RJTrans = BL.*([0.2021 0.1055 -0.1477] );
        LJTrans = BL.*([0.2021 -0.1055 -0.1477]);

        %Assemble full state
        alpha = pp(1);

        p = pfull;
        p(idx) = pp(2:end);

        bodyscrew = p(1:6);
        Gbody = screw2homo(bodyscrew);
        Tbody = Gbody(1:3,4);
        Rbody = Gbody(1:3,1:3);

        %Apply extra roll rotation
        Rbody = Rx(alpha)*Rbody;
        Gbody = [Rx(alpha) zeros(3,1); zeros(1,3) 1]*Gbody;

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

        %===========================================================
        % Determine the wing vectors

        % Vector from wing hinge to wing tip
        Vr = GR(1:3,2);
        Vl = -GL(1:3,2);
        %Make them point in same direction
        Vl(2) = -Vl(2);

        F = norm(Vl - Vr);

        %keyboard


    end
    function [c,ceq] = cfun(pp)
        persistent xorg yorg zorg

        c = [];

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
        alpha = pp(1);

        p = pfull;
        p(idx) = pp(2:end);

        bodyscrew = p(1:6);
        Gbody = screw2homo(bodyscrew);
        Tbody = Gbody(1:3,4);
        Rbody = Gbody(1:3,1:3);

        %Apply extra roll rotation
        Rbody = Rx(alpha)*Rbody;
        
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


        ceq = NewWingPts - WingPts;
        
    end

end