function [camnew,avgres] = dltfu_nonlin(F,L,Cut)
% Description:  Program to calculate DLT coefficient for one camera
%               Note that at least 6 (valid) calibration points are needed
%               function [A,avgres] = dltfu(F,L,Cut)
% Input:        - F      matrix containing the global coordinates (X,Y,Z)
%                        of the calibration frame
%                        e.g.: [0 0 20;0 0 50;0 0 100;0 60 20 ...]
%               - L      matrix containing 2d coordinates of calibration 
%                        points seen in camera (same sequence as in F)
%                        e.g.: [1200 1040; 1200 1360; ...]
%               - Cut    points that are not visible in camera;
%                        not being used to calculate DLT coefficient
%                        e.g.: [1 7] -> calibration point 1 and 7 
%                        will be discarded.
%                        This input is optional (default Cut=[]) 
% Output:       - A      11 DLT coefficients
%               - avgres average residuals (measure for fit of dlt)
%                        given in units of camera coordinates
%
% Author:       Christoph Reinschmidt, HPL, The University of Calgary
% Date:         January, 1994
% Last changes: November 29, 1996
% Version:      1.0
% References:   Woltring and Huiskes (1990) Stereophotogrammetry. In
%               Biomechanics of Human Movement (Edited by Berme and
%               Cappozzo). pp. 108-127.




if nargin==2; Cut=[]; end;

if size(F,1) ~= size(L,1)
disp('# of calibration points entered and seen in camera do not agree'), return
end

m=size(F,1); Lt=L'; C=Lt(:);

for i=1:m
  B(2*i-1,1)  = F(i,1); 
  B(2*i-1,2)  = F(i,2); 
  B(2*i-1,3)  = F(i,3);
  B(2*i-1,4)  = 1;
  B(2*i-1,9)  =-F(i,1)*L(i,1);
  B(2*i-1,10) =-F(i,2)*L(i,1);
  B(2*i-1,11) =-F(i,3)*L(i,1);
  B(2*i,5)    = F(i,1);
  B(2*i,6)    = F(i,2);
  B(2*i,7)    = F(i,3);
  B(2*i,8)    = 1;
  B(2*i,9)  =-F(i,1)*L(i,2);
  B(2*i,10) =-F(i,2)*L(i,2);
  B(2*i,11) =-F(i,3)*L(i,2);
end

% Cut the lines out of B and C including the control points to be discarded
Cutlines=[Cut.*2-1, Cut.*2];
B([Cutlines],:)=[];
C([Cutlines],:)=[];

% Solution for the coefficients
A=pinv(B)*C;

% Initial guess
A0 = A(2:11);

%Image center
imgres = [512 512];
u0 = (imgres(1) - 1)/2;
v0 = (imgres(2) - 1)/2;

% options = optimset('display','iter','gradobj','off','GradConstr','on');
%Anew = fmincon(@errfun,A0,[],[],[],[],A0-delta,A0+delta,@constfun,options);
options = optimset('maxfunevals',1e9,'MaxIter',1e6);
%Anew = fminunc(@errfun,A0,options);
Anew = lsqnonlin(@errfun1,A0,[],[],options);
Anew = hh1(Anew);

%Calculate the extrinsic parameters from
cam = dlt2cam(Anew);


%Let's transform the world coordinates into the camera frame
campts = [cam.R -cam.R*cam.T]*[F ones(size(F,1),1)]';
%Normalize the point coordinates by perspective.
xn = campts(1:2,:)./[campts(3,:);campts(3,:)];
xn = xn';

%Let's get initial estimate for u-v-scale
A1 = A0;
Anew = lsqnonlin(@errfun2,A1,[],[],options);
Anew = hh2(Anew);

%Calculate the extrinsic parameters from
cam1 = dlt2cam(Anew);

%Now minimize the error to estimate the internal parameters of the camera
%P = [du dv u0 v0 k1];
P0 = [-cam.f cam.u0 0];

options = optimset('maxfunevals',1e9,'MaxIter',1e6,'display','iter');
[Pnew,avgres] = lsqnonlin(@internalerrfun,P0,[],[],options);

camnew.T = cam.T;
camnew.R = cam.R;
camnew.u0 = Pnew(3:4);
camnew.f = Pnew(1:2);
camnew.k = Pnew(5);

% %Plot Reprojected points
% % Solution for the coefficients
% D=B*Anew;
c = reshape(C,2,[]);
% d = reshape(D,2,[]);
d = projectpts(F,camnew)';
figure; ploth = plot(c(1,:),c(2,:),'*',d(1,:),d(2,:),'.'); axis equal
for i = 1:size(c,2)
    hold on;
    text(c(1,i),c(2,i),num2str(i));
end
legend(ploth,'original','estimated');


%keyboard
    function [f] = errfun1(A)
        %Incoporate nonlinear constraint
        a = hh1(A);
        %f = 0.5*(sum( (B*a - C).^2));
        f = (B*a - C);
    end

    function [f] = errfun2(A)
        %Incoporate nonlinear constraint
        a = hh2(A);
        %f = 0.5*(sum( (B*a - C).^2));
        f = (B*a - C);
    end

    function [f] = internalerrfun(P)
        r2 = sum(xn.^2,2);
        
        k1 = P(5);
        KK = [P(1) 0 P(3)
            0 P(2) P(4)
            0 0 1];
        xd = xn.* repmat(1 + k1*r2,1,2);
        
        xbar = KK*[xd ones(size(xd,1),1)]';
        
        Cmod = reshape(xbar(1:2,:),[],1);
        
        %f = 0.5*(sum( (B*a - C).^2));
        f = (Cmod - C);
    end

    function a = hh1(A)
        A = [0;A];
        
        %This is the normal constraint calculated from the Modified DLT approach 
        a(1,:) = (-(A(11)*A(2)-A(10)*A(3))*(A(11)*A(6)-A(10)*A(7))...
           + (A(10)*A(2)+A(11)*A(3))*A(5)*A(9) - (A(2)*A(6)+A(3)*A(7))*A(9)^2)/...
           ((A(10)^2 + A(11)^2)*A(5) - (A(10)*A(6)+A(11)*A(7))*A(9));
                
        a(2:11,:) = A(2:11,:);
    end

    function a = hh2(A)
        A = [0;A];
        
        D2 = 1 / (A(10)^2+A(11)^2+A(9)^2);

        %This constraint fixes the principal point at the center of the
        %image axis.
        a(1,:) = (u0*v0 - D2*(A(2)*A(6)+A(3)*A(7))) / (D2*A(5));
        
        a(2:11,:) = A(2:11,:);
    end

end


    
    