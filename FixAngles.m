% Fix alpha Values
obj = {'body','left_wing','right_wing'};

load kine/SavedKinematics/Wings/exp101_ebraheem.mat
for i = 1:length(obj)
    cur_obj = obj{i};
    for frame = 1:258
        
        alpha = (pi/180) * data.kine.(cur_obj).data.params(1,frame);%Get stored parameter value       
        if alpha < 0
            alpha = 2*pi - alpha;
        end
        c1 = data.kine.(cur_obj).data.coords(1,:,frame);
        c2 = data.kine.(cur_obj).data.coords(2,:,frame);
        if i == 1
            dc = c1 - c2;
        else
            dc = c2 - c1;
        end
        
        [psi,theta,c_length] = cart2sph(dc(1),dc(2),dc(3));
        
        alpha = alpha;%pi-alpha;%(2*pi)-alpha; %DEBUG NEW SYSTEM: for new will just be alpha
        theta = -theta; % reverse theta because cart2sph defines elevation above the xy-plane as positive, but this is a negative rotation around the y-axis
        psi = psi; % don't reverse because rotate from body to lab frame

        % Calculate quaternion, rotate, scale, translate
        q = eulzyx2quat(alpha,theta,psi);
        
        data.kine.(cur_obj).data.eulzyx(:,frame) = [alpha;theta;psi];
        data.kine.(cur_obj).data.quat(:,frame) = q;
    end
end
datanew = data;
save('kine/SavedKinematics/Wings/exp101_ebraheemNEW.mat','datanew');
