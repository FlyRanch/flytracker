% view Image from 3 camera views.

% This is a first attempt to make a 3D model from the data given to us by
% Will Dickson

%Load Camera Data
load DLT_coeff

% The joint of the fly wing is located at [0 -1.3412 0.1442]
XTrans = [0 -1.3412 0.1442];
load flywing
load flybody

body = flybody;
wing = flywing;
wing.pts = wing.pts - repmat(XTrans,size(wing.pts,1),1);

%Approximate translation from center of fly body to the wing joint location
%JTrans = [0 +/-0.35 0.5]
RJTrans = [0 -0.36 0.5];
LJTrans = [0 0.36 0.5];

bodyscale = 1;
wingscale = 1;

%Body
body.pts = bodyscale.*body.pts;

%Left Wing
Lwing = wing;
Lwing.pts = wingscale.*Lwing.pts;
Lwing.pts = Lwing.pts + repmat(bodyscale.*LJTrans,size(Lwing.pts,1),1);

%Right Wing
%Rotate about z-axis by 180 degrees
Rz = @(theta) [cos(theta) -sin(theta) 0
                sin(theta) cos(theta) 0
                0 0 1];
Rwing = wing;
Rwing.pts = wingscale.*Rwing.pts;
Rwing.pts = (Rz(pi)*Rwing.pts')' + repmat(bodyscale.*RJTrans,size(Rwing.pts,1),1);

figure; 
trisurf(body.TRI,body.pts(:,1),body.pts(:,2),body.pts(:,3),'facecolor','green');
hold on;
trisurf(Rwing.TRI,Rwing.pts(:,1),Rwing.pts(:,2),Rwing.pts(:,3),'facecolor','blue');
trisurf(Lwing.TRI,Lwing.pts(:,1),Lwing.pts(:,2),Lwing.pts(:,3),'facecolor','red');
axis equal

%get image points
figure;
for i = 1:3
    rw(i).pts = dlt_3D_to_2D(DLT(i).coeff,Rwing.pts);
    lw(i).pts = dlt_3D_to_2D(DLT(i).coeff,Lwing.pts);
    bd(i).pts = dlt_3D_to_2D(DLT(i).coeff,body.pts);
    
    subplot(1,3,i), hold on, plot(rw(i).pts(:,1),rw(i).pts(:,2),'b.');
    subplot(1,3,i), hold on, plot(lw(i).pts(:,1),lw(i).pts(:,2),'r.');
    subplot(1,3,i), hold on, plot(bd(i).pts(:,1),bd(i).pts(:,2),'g.');
end
