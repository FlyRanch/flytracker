function [Gbody,GL,GR] = GTransform(p,PAR)

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