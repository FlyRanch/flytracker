function calc_angles_stroke_plane_centered(frame)


global KineData

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

objects = get(data.handles.cur_object_menu,'String');

% BEGIN CONSTANT DECLARATIONS ==========================
% Kine_v1_0 constants
head=1;tail=2;rwh=3;rwt=4;lwh=5;lwt=6;rwn=7;lwn=8;
right=1;left=2;

% Kine_v2_1 constants
left_wing=objects{1};body=objects{2};right_wing=objects{3};
head=1;tail=2;
hinge=1;tip=2;
axis_name = ['x' 'y' 'z'];
phi = 1; theta = 2;
% END CONSTANT DECLARATIONS ============================


% Initialize 'KineData' struct
KineData.pos.x(:,frame) = zeros(8,1);
KineData.pos.y(:,frame) = zeros(8,1);
KineData.pos.z(:,frame) = zeros(8,1);

% THE FOLLOWING IS TEMPORARY CODE!!  UPDATE ASAP
%Hovering
% KineData.StrokePlaneDeg = -50;
%Left
KineData.StrokePlaneDeg = -60;
% KineData.StrokePlaneDeg = -54;
% END TEMPORARY CODE

% Transfer data in 'data' struct to 'KineData' struct
for axis = 1:3
	KineData.pos.(axis_name(axis))(head, frame) = data.kine.(body).data.coords(head,axis,frame);
	KineData.pos.(axis_name(axis))(tail, frame) = data.kine.(body).data.coords(tail,axis,frame);
	KineData.pos.(axis_name(axis))(rwh, frame) = data.kine.(right_wing).data.coords(hinge,axis,frame);
	KineData.pos.(axis_name(axis))(rwt, frame) = data.kine.(right_wing).data.coords(tip,axis,frame);
	KineData.pos.(axis_name(axis))(lwh, frame) = data.kine.(left_wing).data.coords(hinge,axis,frame);
	KineData.pos.(axis_name(axis))(lwt, frame) = data.kine.(left_wing).data.coords(tip,axis,frame);
end

KineData.pos.angles(right,frame) = deg2rad(data.kine.(right_wing).data.params(1, frame));
KineData.pos.angles(left,frame) = deg2rad(data.kine.(left_wing).data.params(1, frame));

% Execute Kine_v1_0 code
calcBodyStrokePlaneCentred(frame)

% Extract calculated data back into 'data' struct
for axis = 1:3
	data.kine.(body).BCdata.coords(head,axis,frame) = KineData.posBC.(axis_name(axis))(head,frame);
	data.kine.(body).BCdata.coords(tail,axis,frame) = KineData.posBC.(axis_name(axis))(tail,frame);
	data.kine.(right_wing).BCdata.coords(hinge,axis,frame) = KineData.posBC.(axis_name(axis))(rwh,frame);
	data.kine.(right_wing).BCdata.coords(tip,axis,frame) = KineData.posBC.(axis_name(axis))(rwt,frame);
	data.kine.(left_wing).BCdata.coords(hinge,axis,frame) = KineData.posBC.(axis_name(axis))(lwh,frame);
	data.kine.(left_wing).BCdata.coords(tip,axis,frame) = KineData.posBC.(axis_name(axis))(lwt,frame);
end

data.kine.(left_wing).BCdata.params(1,frame) = KineData.robo.alpha.left(frame);
data.kine.(right_wing).BCdata.params(1,frame) = KineData.robo.alpha.right(frame);

data.kine.(left_wing).BCdata.angles(1,theta,frame) = KineData.robo.theta.left(frame);
data.kine.(left_wing).BCdata.angles(1,phi,frame) = KineData.robo.phi.left(frame);

data.kine.(right_wing).BCdata.angles(1,theta,frame) = KineData.robo.theta.right(frame);
data.kine.(right_wing).BCdata.angles(1,phi,frame) = KineData.robo.phi.right(frame);

% % Calculate angles for body centered coordinates
% coord_diff = data.kine.(left_wing).BCdata.coords(2,:,frame) - data.kine.(left_wing).BCdata.coords(1,:,frame);
% [theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
% data.kine.(left_wing).BCdata.angles(:, :, frame) = [phi theta r];
% 
% coord_diff = data.kine.(right_wing).BCdata.coords(2,:,frame) - data.kine.(right_wing).BCdata.coords(1,:,frame);
% [theta phi r] = cart2sph(coord_diff(1), coord_diff(2), coord_diff(3));
% data.kine.(right_wing).BCdata.angles(:, :, frame) = [phi theta r];

% Save 'data' struct
guidata(controller,data);