function xi = homo2screw(G)

% Convert an associated homogeneous transformation to a screw motion
%
% xi = [v ; w]
%
% G = [R T
%      0 1];
%
% Code also taken from 'rodrigues.m' written by Pietro Perona and used in
% Caltech camera calibration code

R = G(1:3,1:3);
T = G(1:3,4);


%==================================
% First use routine to determine omega axis from rotation matrix
%
% project the rotation matrix to SO(3);
[U,S,V] = svd(R);
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


    omega = om*theta;

    domegadvar2 = [theta*eye(3) om];
    dout = domegadvar2 * dvar2dvar * dvardR;

else
    if tr > 0; 			% case norm(om)=0;

        omega = [0 0 0]';

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

            omega = theta * [uabs; vabs; wabs] .* svec;

        end;
    end;
end;

% Next, use Rodrigues from the calculated omega 
theta = norm(omega);

%make omega unit length
omega = omega./theta;

if theta == 0;
    xi(1:3,:) = T;
    xi(4:6,:) = zeros(3,1);
else
    beta = sin(theta);
    gamma = 1-cos(theta);
    omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        
    R = eye(3) + omegav*beta + omegav*omegav*gamma;
    
    
    A = (eye(3) - R)*omegav + omega*omega'.*theta;
    
    v = inv(A)*T;
    
    xi(1:3,:) = v(:);
    xi(4:6,:) = omega(:)*theta;
    
end