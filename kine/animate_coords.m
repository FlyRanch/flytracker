% Script to animate body coords in 3D

function animate_coords(exp)

switch exp
    
    case 0 
        
        disp('Options: 79 | 82 | 83 | 86 | 88')
        return
        
    case 79
        load bodykin_exp079
        frames = 231:372;
        
    case 82
        load bodykin_exp082
        frames = 1:162;
        
    case 83
        load bodykin_exp083
        frames = 367:588;
        
    case 86
        load bodykin_exp086
        frames = 264:464;
        
    case 88
        load bodykin_exp088
        frames = 9:130;
        
    case 89
        load bodykin_exp089
        frames = 42:117;
        
    case 103
        load bodykin_exp103
        frames = 233:442;

    otherwise
        warndlg('No Data for that experiment yet')
        return
        
end


% Get axes ranges
max_vals = max(data.kine.body.coords(:,:,frames),[],1);
max_vals = max(max_vals,[],3);

min_vals = min(data.kine.body.coords(:,:,frames),[],1);
min_vals = min(min_vals,[],3);

marg = 3;
max_vals = max_vals + marg;
min_vals = min_vals - marg;

figure
ax = axes;
hold on
set(ax,'NextPlot','add')
view(-130,30)
% view(-172,4)
set(ax,'XLim',[min_vals(1) max_vals(1)],'YLim',[min_vals(2) max_vals(2)],'ZLim',[min_vals(3) max_vals(3)])
set(ax,'ZDir','reverse')
set(ax,'XGrid','on','YGrid','on','ZGrid','on')

title('Press any key to start the animation')

for frame = frames
    
    x = data.kine.body.coords(:,1,frame);
    y = data.kine.body.coords(:,2,frame);
    z = data.kine.body.coords(:,3,frame);
    
    head_x = data.kine.body.coords(1,1,frame);
    head_y = data.kine.body.coords(1,2,frame);
    head_z = data.kine.body.coords(1,3,frame);
    
    old = findobj('Type','line');
    set(old,'Color','c')
    
    plot3(x, y, z, 'b-', 'LineWidth', 2, 'EraseMode','none');
    plot3(head_x, head_y, head_z, 'bo','EraseMode','none');
    xlabel(['Frame ',num2str(frame)])
    
    if frame == frames(1)
        pause
    end
    
    pause(0.02)
    
end



% %--------------------------------------------------------------------------
% function data_out = smooth_data(data_in)
% 
% x = data_in(1,:,:);
% y = data_in(2,:,:);
% z = data_in(3,:,:);
% 
% f = 1:size(data_in,[],3);
% 
% xx = spline(f,x,