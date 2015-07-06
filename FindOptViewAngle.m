function cam = FindOptViewAngle(a,L,F,cam,PAR)

for i = 1:PAR.numcam
    figure(i);
    clf;
    
%     plot(L(:,1,i),L(:,2,i),'g*');
%     ax1 = gca;
    

    plot3(F(:,1),F(:,2),F(:,3),'k.');
    axis equal
    ax(i) = gca;
    
    %Determine the point in space that the camera is looking at.
    cam(i).xyz = dlt_2D_to_3Diter(a(:,i),PAR.imgres(1)/2,PAR.imgres(2)/2);
    
    %Set the camera parameters to get correct perspective transform
    set(ax(i),'CameraPosition',cam(i).C,'CameraUpVector',cam(i).R(2,:),...
        'Projection','perspective','CameraTarget',cam(i).xyz);
            
    ang0 = get(ax(i),'CameraViewAngle');
    
    
 
    options = optimset('display','none');
    %cam(i).newxyz = lsqnonlin(@viewerr,cam(i).xyz,[],[],options);
    %cam(i).viewangle = lsqnonlin(@viewerr,ang0,[],[],options);
    cam(i).viewangle = fminbnd(@viewerr,ang0,ang0+10,options);
    
%     figure(i+20);
%     clf;
%     plot(L(:,1,i),L(:,2,i),'*'); 
%     axis(0.5+[0 PAR.imgres(1) 0 PAR.imgres(2)]);
%     
%     hold on; plot(npts(:,1),npts(:,2),'r.')
%     hold on; plot(PAR.imgres(1)/2,PAR.imgres(2)/2,'ko')
%     keyboard
    
end
close all

    function F = viewerr(ang)
        set(ax(i),'CameraPosition',cam(i).C,'CameraUpVector',cam(i).R(2,:),...
            'Projection','perspective','CameraTarget',cam(i).xyz,'CameraViewAngle',ang);
        
        imnew = rendercalib(i,PAR.imgres);
        [r,c] = find(imnew == 255);
        
        npts = [c PAR.imgres(2)-r];
        ii = kdtreeidx_nrml(npts,L(:,:,i));
        
        F = npts(ii,:) - L(:,:,i);
        
        F = reshape(F',[],1);
        
        F = sum(F.^2);
    end
    
end
