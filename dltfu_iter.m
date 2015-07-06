function [Anew,avgres] = dltfu_iter(F,L,Cut)
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
D=B*A;
R=C-D;
avgres=0.5*(sum( (B*A - C).^2));

delta = [5*ones(3,1);0;5*ones(3,1);0;5*ones(3,1)];
% Initial guess
A0 = A(2:11);

%Image center
imgres = [512 512];
u0 = (imgres(1) - 1)/2;
v0 = (imgres(2) - 1)/2;

% options = optimset('display','iter','gradobj','off','GradConstr','on');
%Anew = fmincon(@errfun,A0,[],[],[],[],A0-delta,A0+delta,@constfun,options);
options = optimset('maxfunevals',1e9,'MaxIter',1e6,'display','none');
%Anew = fminunc(@errfun,A0,options);
[Anew,avgres] = lsqnonlin(@errfun,A0,[],[],options);
Anew = hh(Anew);

%Plot Reprojected points
% Solution for the coefficients
D=B*Anew;
c = reshape(C,2,[]);
d = reshape(D,2,[]);
figure; ploth = plot(c(1,:),c(2,:),'*',d(1,:),d(2,:),'.'); axis equal
for i = 1:size(c,2)
    hold on;
    text(c(1,i),c(2,i),num2str(i));
end
legend(ploth,'original','estimated');


%keyboard
    function [f] = errfun(A)
        %Incoporate nonlinear constraint
        a = hh(A);
        %f = 0.5*(sum( (B*a - C).^2));
        f = (B*a - C);
    end

    function a = hh(A)
        A = [0;A];
        
        %This is the normal constraint calculated from the Modified DLT approach 
        a(1,:) = (-(A(11)*A(2)-A(10)*A(3))*(A(11)*A(6)-A(10)*A(7))...
           + (A(10)*A(2)+A(11)*A(3))*A(5)*A(9) - (A(2)*A(6)+A(3)*A(7))*A(9)^2)/...
           ((A(10)^2 + A(11)^2)*A(5) - (A(10)*A(6)+A(11)*A(7))*A(9));
        
        %This constraint fixes the principal point at the center of the
        %image axis.
%         a(1,:) = A(9)*(u0 + (A(11)*A(6) - A(10)*A(7))*(-A(3)+A(11)*u0) / ...
%           (A(10)^2*A(6) + A(10)*A(11)*A(7) + A(6)*A(9)^2 - A(10)*...
%           (A(10)^2+A(11)^2+A(9)^2)*v0) );
        
        a(2:11,:) = A(2:11,:);
    end


    function [c,ceq,gc,gceq] = constfun(A);
        c = [];
        
        ceq = (A(1)*A(9) + A(2)*A(10) + A(3)*A(11)) * ...
            (A(5)*A(9) + A(6)*A(10) + A(7)*A(11)) / (A(9)^2 + A(10)^2 + A(11)^2) ...
            - (A(1)*A(5) + A(2)*A(6) + A(3)*A(7));
        
        if nargout > 2
            %These are calculated from the mathematica notebook
            %ModifiedDLT.nb
            
            gc = [];
            
            t1 = A(1)*A(9) + A(2)*A(10) + A(3)*A(11);
            t2 = A(5)*A(9) + A(6)*A(10) + A(7)*A(11);
            t3 = A(9)^2 + A(10)^2 + A(11)^2;
            
            gceq(1,:) = -A(5) + A(9)*t2 / t3;
            
            gceq(2,:) = -A(6) + A(10)*t2 / t3;
            
            gceq(3,:) = -A(7) + A(11)*t2 / t3;
            
            gceq(4,:) = 0;
            
            gceq(5,:) = -A(1) + A(9)*t1 / t3;
            
            gceq(6,:) = -A(2) + A(10)*t1 / t3;
            
            gceq(7,:) = -A(3) + A(11)*t1 / t3;
            
            gceq(8,:) = 0;
            
            gceq(9,:) = ( -2*A(9)*t1*t2 + A(5)*t1*t3 + A(1)*t2*t3) / t3^2;
            
            gceq(10,:) = ( -2*A(10)*t1*t2 + A(6)*t1*t3 + A(2)*t2*t3) / t3^2;
            
            gceq(11,:) = ( -2*A(11)*t1*t2 + A(7)*t1*t3 + A(3)*t2*t3) / t3^2;
            
        end
            
    end
end


    
    