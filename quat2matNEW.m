function Y = quat2matNEW(X)
% Y = quat2mat(X) converts rotations in SO(3) between matrix and
% quaternion representation.  If one in input, the other is
% output.
% Use Rodriguez formula to convert to angle axis, and then angle/axis to
% quaterion

eps = realmin;
dim = prod(size(X));
switch dim
    case 9 %A rotation matrix is the input
        % project the rotation matrix to SO(3);
        [U,S,V] = svd(X);
        R = U*V';

        tr = (trace(R)-1)/2;
        dtrdR = [1 0 0 0 1 0 0 0 1]/2;
        theta = real(acos(tr));


        if sin(theta) >= 1e-4,

            dthetadtr = -1/sqrt(1-tr^2);

            dthetadR = dthetadtr * dtrdR;
            % var1 = [vth;theta];
            vth = 1/(2*sin(theta));
            dvthdtheta = -vth*cos(theta)/sin(theta);
            dvar1dtheta = [dvthdtheta;1];

            dvar1dR =  dvar1dtheta * dthetadR;


            om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';

            dom1dR = [0 0 0 0 0 1 0 -1 0;
                0 0 -1 0 0 0 1 0 0;
                0 1 0 -1 0 0 0 0 0];

            % var = [om1;vth;theta];
            dvardR = [dom1dR;dvar1dR];

            % var2 = [om;theta];
            om = vth*om1;
            domdvar = [vth*eye(3) om1 zeros(3,1)];
            dthetadvar = [0 0 0 0 1];
            dvar2dvar = [domdvar;dthetadvar];


            out = om*theta;

        else
            if tr > 0; 			% case norm(om)=0;

                out = [0 0 0]';

                dout = [0 0 0 0 0 1/2 0 -1/2 0;
                    0 0 -1/2 0 0 0 1/2 0 0;
                    0 1/2 0 -1/2 0 0 0 0 0];
            else

                % case norm(om)=pi;
                if(0)

                    %% fixed April 6th by Bouguet -- not working in all cases!
                    out = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
                    %keyboard;

                else

                    % Solution by Mike Burl on Feb 27, 2007
                    % This is a better way to determine the signs of the
                    % entries of the rotation vector using a hash table on all
                    % the combinations of signs of a pairs of products (in the
                    % rotation matrix)

                    % Define hashvec and Smat
                    hashvec = [0; -1; -3; -9; 9; 3; 1; 13; 5; -7; -11];
                    Smat = [1,1,1; 1,0,-1; 0,1,-1; 1,-1,0; 1,1,0; 0,1,1; 1,0,1; 1,1,1; 1,1,-1;
                        1,-1,-1; 1,-1,1];

                    M = (R+eye(3,3))/2;
                    uabs = sqrt(M(1,1));
                    vabs = sqrt(M(2,2));
                    wabs = sqrt(M(3,3));

                    mvec = [M(1,2), M(2,3), M(1,3)];
                    syn  = ((mvec > 1e-4) - (mvec < -1e-4)); % robust sign() function
                    hash = syn * [9; 3; 1];
                    idx = find(hash == hashvec);
                    svec = Smat(idx,:)';

                    out = theta * [uabs; vabs; wabs] .* svec;

                end                
            end
        end
        out = out(:);
        theta = norm(out);
        omega = out./theta;
        Y = [omega.*sin(theta/2) ; cos(theta/2)];

    case 4 %A quaternion is the input
        % Check to make sure that it is unit magnitude
        %   if norm(X) ~= 1
        %     error(['the quaternion does not have unit length.  Length = ' ...
        % 	   'num2str(norm(X))']);
        %   end
        
        e_0 = X(4);
        e_x = X(1);
        e_y = X(2);
        e_z = X(3);

        Y = [ e_0^2+e_x^2-e_y^2-e_z^2 2*e_x*e_y-2*e_0*e_z     2*e_x*e_z+2*e_0*e_y;...
            2*e_x*e_y+2*e_0*e_z     e_0^2-e_x^2+e_y^2-e_z^2 2*e_y*e_z-2*e_0*e_x;...
            2*e_x*e_z-2*e_0*e_y     2*e_y*e_z+2*e_0*e_x     e_0^2-e_x^2-e_y^2+e_z^2 ];

%         w = X(4);
%         x = X(1);
%         y = X(2);
%         z = X(3);
% 
%         Y = [w^2 + x^2 - y^2 - z^2, 2*x*y+2*w*z, 2*x*z-2*w*y
%             2*x*y-2*w*z, w^2 - x^2 + y^2 - z^2, 2*y*z+2*w*x
%             2*x*z+2*w*y, 2*y*z-2*w*x, w^2 - x^2 - y^2 + z^2];
        
end


